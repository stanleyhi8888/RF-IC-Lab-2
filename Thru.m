
% File I/O
thru=sparameters('q2_de.s2p');
L1meas=sparameters('q2.s2p');

z0=50;
freq=thru.Frequencies;

%Data transform
thru_abcd=s2abcd(thru.Parameters);%get Zs Yp
for i=1:size(freq)
    Tmeas(:,:,i)=s2abcd(L1meas.Parameters(:,:,i));
    Zs(i)=thru_abcd(1,2,i)/2;
    Yp(i)=(thru_abcd(1,1,i)-1)/(thru_abcd(1,2,i));
end

%acquire TL TR
for x=1:size(freq)
    for y=1:2
        for z=1:2
            if z==1&&y==1
            TL(z,y,x)=1;
            TR(z,y,x)=1+Yp(1,x)*Zs(1,x);
            elseif (z==1)&&(y==2)
            TL(z,y,x)=Zs(1,x);
            TR(z,y,x)=Zs(1,x);
            elseif z==2&&y==1
            TL(z,y,x)=Yp(1,x);
            TR(z,y,x)=Yp(1,x);
            elseif z==2&&y==2
            TL(z,y,x)=1+Yp(1,x)*Zs(1,x);
            TR(z,y,x)=1;
            else
            print("e04");
            break
            end
        end
    end
end

%get SDUT
for i=1:size(freq)
    TL_inv(:,:,i)=inv(TL(:,:,i));
    TR_inv(:,:,i)=inv(TR(:,:,i));
    TDUT(:,:,i)=TL_inv(:,:,i)*Tmeas(:,:,i)*TR_inv(:,:,i);
    SDUT(:,:,i)=abcd2s(TDUT(:,:,i));
    YDUT(:,:,i)=abcd2y(TDUT(:,:,i));
    %TLine
    E(i)=180*acosh(TDUT(1,1,i))/(pi);
    Z(i)=z0*sqrt(((1+SDUT(1,1,i))^2-SDUT(2,1,i)^2)/((1-SDUT(1,1,i))^2-SDUT(2,1,i)^2));
    %Resistor
    %R(i)=1/real(YDUT(1,1,i));
    %Im(i)=imag(YDUT(1,1,i)+YDUT(2,2,i))/2;
end

rfwrite(SDUT(:,:,1:length(freq)),freq,'DUT.s2p');
DUT=sparameters('DUT.s2p');

%{
TRTL_cal=zeros(2,2,length(freq));
%verification
for i=1:length(freq) 
    TRTL_cal(1:2,1:2,i)=[1+2*Yp(i)*Zs(i) 2*Zs(i); 2*Yp(i)*(Yp(i)*Zs(i)+1) 1+2*Yp(i)*Zs(i)];
end
%}

%-----------------fitting process (禁術)--------------------------
%for TLine
%{
exp_coeff=0.3;
scal_coeff=21;
for i=1:size(freq)
    Z(i)=scal_coeff*Z(i)^(exp_coeff);
end
%}
%for Resistor
%{
exp_coeff=0.15;
scal_coeff=26;
for i=1:size(freq)
    R(i)=scal_coeff*R(i)^(exp_coeff);
end
%}


start=51;
stop=100;

%q3 varient
figure
plot(freq(start:stop),Z(start:stop));
figure
plot(freq(start:stop),E(start:stop));

%q4 varient
%{
figure
plot(freq(start:stop),R(start:stop));
figure
plot(freq(start:stop),Im(start:stop));
%}
