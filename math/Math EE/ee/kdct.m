% Usage: Y = kdct(y).  Substitute for Matlab's dct() command. Not as 
% efficient, but it works. If input y is row vector, transforms the row, 
% otherwise transforms column by column.
function Y = kdct(y)

[M,N] = size(y);

if M==1 %Do a row vector transform
    y2 = [y, y(:,N:-1:1)]; %Extend vector to double length
    Y2 = fft(double(y2)); %Let Matlab do the dirty work
    Y = real(Y2(1:N).*exp(-0.5*pi*1i*(0:N-1)/N))/sqrt(2*N); %Restrict, weight
    Y(1) = Y(1)/sqrt(2);
    
else %Do a matrix (or column vector) column by column
    y2 = [y; y(M:-1:1,:)];  %Reflect each column
    Y2 = fft(double(y2)); %Works on each column
    Y = zeros(M,N);
    for k=1:M
        r = exp(-0.5*pi*1i*(k-1)/M)/sqrt(2*M);
        Y(k,:) = real(Y2(k,:)*r);
    end
    Y(1,:) = Y(1,:)/sqrt(2);
end
