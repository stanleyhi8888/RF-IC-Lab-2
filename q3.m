%File I/O
SL3meas=sparameters('q3_5.8cm.s2p');
SL1meas=sparameters('q2.s2p');

z0=50;
freq=SL1meas.Frequencies;

%StoT TRANSFORM
for i=1:length(freq)
    TL1meas_inv(:,:,i)=inv(s2abcd(SL1meas.Parameters(:,:,i)));
    TL3meas(:,:,i)=s2abcd(SL3meas.Parameters(:,:,i));
    TL_mul_TR(:,:,i)= TL3meas(:,:,i)*TL1meas_inv(:,:,i)*TL3meas(:,:,i);
    Zs(i)=TL_mul_TR(1,2,i)/2;
    Yp(i)=(TL_mul_TR(1,1,i)-1)/(TL_mul_TR(1,2,i));
end

for x=1:length(freq)
    for y=1:2
        for z=1:2
            if z==1&&y==1
            TL(z,y,x)=1;
            TR(z,y,x)=1+Yp(i)*Zs(i);
            elseif (z==1)&&(y==2)
            TL(z,y,x)=Zs(i);
            TR(z,y,x)=Zs(i);
            elseif z==2&&y==1
            TL(z,y,x)=Yp(i);
            TR(z,y,x)=Yp(i);
            elseif z==2&&y==2
            TL(z,y,x)=1+Yp(i)*Zs(i);
            TR(z,y,x)=1;
            else
            print("e04");
            break
            end
        end
    end
end

%get SDUT
for i=1:length(freq)
    TL_inv(:,:,i)=inv(TL(:,:,i));
    TR_inv(:,:,i)=inv(TR(:,:,i));
    TDUT(:,:,i)=TL_inv(:,:,i)*Tmeas(:,:,i)*TR_inv(:,:,i);
    SDUT(:,:,i)=abcd2s(TDUT(:,:,i));
    %TLine
    E(i)=180*acosh(TDUT(1,1,i))/(pi);
    Z(i)=z0*sqrt(((1+SDUT(1,1,i))^2-SDUT(2,1,i)^2)/((1-SDUT(1,1,i))^2-SDUT(2,1,i)^2));
end
start=51;
stop=100;
figure
plot(freq(start:stop),Z(start:stop));
figure
plot(freq(start:stop),E(start:stop));

rfwrite(SDUT(:,:,1:length(freq)), freq, 'q3dut.s2p');