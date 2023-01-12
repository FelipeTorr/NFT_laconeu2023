classdef Plasticity
	properties
		HGamma
		H
		ratesA
		ratesB
		omega
		Phifft
		N
		nfft
		tp
		gamma
		range
        plasticityTime
        samplingFrequency
	end
	methods
		function obj=Plasticity(config)
			%Values and calculations from 'Assadzadeh, 2018'
			%obj.ratesA=[-1.01;-1.01;-1.25;-1;-0.4;-1.25;-1.25;-1.25];
			obj.ratesA=config.ratesA;
			%obj.ratesB=[10,1.6,0.8,2.7,0.7,1.1,0.55,0.45];
			obj.ratesB=config.ratesB;
			obj.tp=config.tp;
			obj.omega=config.omega*config.samplingFrequency;
			obj.N=config.Nx*config.Ny;
			obj.Phifft=zeros(5,config.samplingFrequency*2);
			obj.samplingFrequency=config.samplingFrequency;
            obj.nfft=config.samplingFrequency;
			obj.gamma=config.gamma;
			obj.range=config.range;
			obj.plasticityTime=linspace(-0.1,0.1,config.samplingFrequency*0.2);
			fprintf('Plasticity: Instatiation completed\n')
		end

		function obj=init(obj)
			%Calculate transfer function spectrums
			obj.H=obj.calculateH();%Verify this
			Gamma=obj.calculateGamma(obj.gamma,obj.range);
			%obj.HGamma=conj(obj.H(1,:)).*conj(Gamma);
		end

		function [mfft,means]=calculatefft(obj,Phi)
			%Populations spectrums
			mfft=zeros(5,obj.nfft);
			means=zeros(5,1);
			for m=0:4
				[Phi_pop,meanPhi]=obj.zScore(Phi(m*obj.N+1:(m+1)*obj.N,:));
				means(m+1)=meanPhi;
				sfft=fft(Phi_pop);
				mfft(m+1,:)=sfft(1:obj.nfft)/(2*obj.nfft);
			end
			
		end

		function [phi1,means]=zscoreGrid(obj,Phi)
			phi1=zeros(5,obj.samplingFrequency*2);
			means=zeros(5,1);
			for m=0:4
				[Phi_pop,meanPhi]=obj.zScore(Phi(m*obj.N+1:(m+1)*obj.N,:));
				means(m+1)=meanPhi;
				phi1(m+1,:)=Phi_pop;
			end
		end
			
		function [y]=xcorrGrid(obj,Phi,h,m,n)
			corrmn=zeros(256,obj.samplingFrequency*2);
			m=m-1;
			n=n-1;
			nxh=numel(h);
			halfnxh=floor(nxh/2);
			nx=obj.samplingFrequency*2;
			for j=1:256
				aux_corrmn=xcorr(Phi(m*obj.N+j,:),h);
				corrmn(j,:)=aux_corrmn(end-nx-halfnxh:end-halfnxh-1).*Phi(n*obj.N+j,:);
			end
			averageCorr=-sum(corrmn,1)/256;
			y=sum(averageCorr)/(obj.samplingFrequency*2);
		end
			
		function [y]=xcovGrid(obj,Phi,h,m,n)
			corrmn=zeros(256,obj.samplingFrequency*2);
			m=m-1;
			n=n-1;
			nxh=numel(h);
			halfnxh=floor(nxh/2);
			nx=obj.samplingFrequency*2;
			for j=1:256
				aux_corrmn=xcorr(Phi(m*obj.N+j,:)-mean(Phi(m*obj.N+j,:)),h);
				corrmn(j,:)=aux_corrmn(end-nx-halfnxh:end-halfnxh-1).*(Phi(n*obj.N+j,:)-mean(Phi(n*obj.N+j,:)));
			end
			averageCorr=-sum(corrmn,1)/256;
			y=sum(averageCorr)/(obj.samplingFrequency*2);
		end

		function [z,meanx]=zScore(obj,x)
			%z-score of 
			meanx=mean(x(:));
			z=(mean(x,1)-meanx)./std(x(:));
		end
		

		function Gamma=calculateGamma(obj,gamma,range)
			k=(2*pi*3e9)./obj.omega;
			Gamma=((1-1j*obj.omega./gamma).^2+k.^2*range^2).^(-1);
		end
		
		function integ=integrand(obj,phia,phib,indexH,flagGamma)
			if flagGamma==1
				integ=real(phia.*conj(phib).*obj.HGamma);
			else
				integ=real(phia.*conj(phib).*conj(obj.H(indexH,:)));
			end
			
		end
		
		function integrall=integrall(obj,integrand)
			integrall=integrand(1)*(obj.omega(2)-obj.omega(1))/(2*pi);
			for i=2:length(integrand)
				integrall=integrall+integrand(i)*(obj.omega(i)-obj.omega(i-1))/(2*pi);
			end
		end
		
		function result=integralTime(obj,phia,phib,h)
			nxh=numel(h);
			halfnxh=floor(nxh/2);
			nx=numel(phib);
			convolution1=xcorr(phia,h);
            result=-sum(convolution1(end-nx-halfnxh:end-halfnxh-1).*phib)/nx;
            %convolution1=conv(phia,h);
            %result=sum(convolution1(halfnxh:halfnxh+nx-1).*phib)/nx;
		end

		function y=integralAverage(obj,x1,x2,h)
			nx1=numel(x1);
			nxh=numel(h);
			y=zeros(1,nxh);
			for m=1:nx1-nxh
				y=y+obj.integralShort(x2(m+1:m+nxh),x1(m+1:m+nxh),h);
			end
			y=sum(y)/(nx1-nxh);
		end

		function y=integralShort(obj,x1,x2,h)
		    nx1=numel(x1);
		    nx2=numel(x2);
		    nxh=numel(h);
		    halfnxh=floor(nxh/2);
		    len_corr=nx1+nx2-1;
		    y=zeros(1,nx1);
		    for n=1:len_corr
		        if n<=nx1-halfnxh
		            for m=1-n:halfnxh
		                y(n)=y(n)+x1(n+m)*h(m+halfnxh)*x2(n);
		            end
		        else
		            for m=-halfnxh+1:halfnxh-n
		                y(n)=y(n)+x1(n+m)*h(m+halfnxh)*x2(n);
		            end
		        end
		    end
		end

		function H=calculateH(obj)
			H=zeros(8,length(obj.plasticityTime));
			H(:,obj.plasticityTime>0)=-obj.ratesB'*exp(-obj.plasticityTime(obj.plasticityTime>0)*1/obj.tp);
			H(:,obj.plasticityTime<=0)=-obj.ratesB'.*(obj.ratesA)*exp(obj.plasticityTime(obj.plasticityTime<=0)*1/obj.tp); 
		end

		function [strengthRatios,means]=update(obj,Phi)
			[phi1,means]=obj.zscoreGrid(Phi);
			strengthRatios=zeros(8,1);
			%s_ee
			strengthRatios(1)=obj.xcorrGrid(Phi,obj.H(1,:),1,1)/(obj.samplingFrequency);
			%s_ei
			strengthRatios(2)=obj.xcorrGrid(Phi,obj.H(2,:),1,2)/(obj.samplingFrequency);
			%s_es
			strengthRatios(3)=obj.xcorrGrid(Phi,obj.H(3,:),1,4)/(obj.samplingFrequency);
			%s_re
			strengthRatios(4)=obj.xcorrGrid(Phi,obj.H(4,:),3,1)/(obj.samplingFrequency);
			%s_rs
			strengthRatios(5)=obj.xcorrGrid(Phi,obj.H(5,:),3,4)/(obj.samplingFrequency);
			%s_se
			strengthRatios(6)=obj.xcorrGrid(Phi,obj.H(6,:),4,1)/(obj.samplingFrequency);
			%s_sr
			strengthRatios(7)=obj.xcorrGrid(Phi,obj.H(7,:),4,3)/(obj.samplingFrequency);
			%s_sn
			strengthRatios(8)=obj.xcorrGrid(Phi,obj.H(8,:),4,5)/(obj.samplingFrequency);
		end

		function [strengthRatios,means]=updateCorr(obj,Phi)
			[phi1,means]=obj.zscoreGrid(Phi);
			strengthRatios=zeros(8,1);
			etat=length(obj.plasticityTime);
			%s_ee
			strengthRatios(1)=obj.integralTime(phi1(1,:),phi1(1,:),obj.H(1,:))/(obj.samplingFrequency);
			%s_ei
			strengthRatios(2)=obj.integralTime(phi1(1,:),phi1(2,:),obj.H(2,:))/(obj.samplingFrequency);
			%s_es
			strengthRatios(3)=obj.integralTime(phi1(1,:),phi1(4,:),obj.H(3,:))/(obj.samplingFrequency);
			%s_re
			strengthRatios(4)=obj.integralTime(phi1(3,:),phi1(1,:),obj.H(4,:))/(obj.samplingFrequency);
			%s_rs
			strengthRatios(5)=obj.integralTime(phi1(3,:),phi1(4,:),obj.H(5,:))/(obj.samplingFrequency);
			%s_se
			strengthRatios(6)=obj.integralTime(phi1(4,:),phi1(1,:),obj.H(6,:))/(obj.samplingFrequency);
			%s_sr
			strengthRatios(7)=obj.integralTime(phi1(4,:),phi1(3,:),obj.H(7,:))/(obj.samplingFrequency);
			%s_sn
			strengthRatios(8)=obj.integralTime(phi1(4,:),phi1(5,:),obj.H(8,:))/(obj.samplingFrequency);
		end

		function [strengthRatios,means]=updateFreq(obj,Phi)
			[mfft,means]=obj.calculatefft(Phi);
			strengthRatios=zeros(8,1);
			%s_ee
			strengthRatios(1)=obj.integrall(obj.integrand(mfft(1,:),mfft(1,:),1,0));
			%s_ei
			strengthRatios(2)=obj.integrall(obj.integrand(mfft(1,:),mfft(2,:),2,0));
			%s_es
			strengthRatios(3)=obj.integrall(obj.integrand(mfft(1,:),mfft(4,:),3,0));
			%s_re
			strengthRatios(4)=obj.integrall(obj.integrand(mfft(3,:),mfft(1,:),4,0));
			%s_rs
			strengthRatios(5)=obj.integrall(obj.integrand(mfft(3,:),mfft(4,:),5,0));
			%s_se
			strengthRatios(6)=obj.integrall(obj.integrand(mfft(4,:),mfft(1,:),6,0));
			%s_sr
			strengthRatios(7)=obj.integrall(obj.integrand(mfft(4,:),mfft(3,:),7,0));
			%s_sn
			strengthRatios(8)=obj.integrall(obj.integrand(mfft(4,:),mfft(5,:),8,0));
		end

		function H=calculateHFreq(obj)
			%STDP spectrum
			H=obj.ratesB'.*(obj.ratesA*obj.tp*(1+1j*obj.omega)+obj.tp*(1-1j*obj.omega))./(1+(obj.omega*obj.tp).^2);
		end
	end
end
