%   The code was written by Jinhui Shi,Kun Zhan
%   $Revision: 1.0.0.0 $  $Date: 2015/01/06 $ 11:37:27 $

%   Reference:
%   K Zhan, J Shi, Q Li, J Teng, M Wang, 
%   "Image segmentation using fast linking SCM," 
%   in Proc. of IJCNN, vol. 25, pp. 2093-2100, 2015.
function Ib = FL_SCM(I) 
    I = GrayStretch(I,0.98);
    Ts = findTh(I)/255;     
    [r,c] = size(I);
    S = im2double(I);
    Sm = max(S(:));
    w = fspecial('gaussian',7,1);
    LAP = abs(imfilter(S,[1 1 1;1 -8 1;1 1 1],'symmetric'));
    bet = 2.3*LAP;
    f = 0.01;
    TTs = Ts/(1 - f);
    dT = 1/255; n = ceil((Sm - TTs)/dT + 1);
    Th = ones(r,c)*Sm;
    Vt = 10000;
    Y = zeros(r,c); YY = Y; U = Y; 
    F = S;    
    for k = 1:n
        L = imfilter(Y,w,'symmetric');
        fire = 1;
        while fire == 1
            Q = Y;
            U = f.*U + F.*(1 + bet.*L);
            Y = double(U > Th);
            if isequal(Q,Y);
                fire = 0;
            else
                L = imfilter(Y,w,'symmetric');
            end
        end
        YY = YY + Y;
        Th = Th - dT + Vt*Y;
    end    
    Ib = YY;
end
function T = findTh(I)   
    I = double(I);
    mi = min(I(:));
    ma = max(I(:));
    L = ma - mi + 1;
    H = zeros(1,L - 2);
    k = 1;
    for i = (mi + 1):(ma - 1)
        B0 = double(I <= i);
        B1 = double(I > i);        
        n0 = sum(B0(:));
        n1 = sum(B1(:));
        I0 = I.*B0;
        I1 = I.*B1;
        u0 = sum(I0(:))/n0;
        u1 = sum(I1(:))/n1;
        s0 = ((I0 - u0).*B0).^2;
        s1 = ((I1 - u1).*B1).^2;
        si0 = sum(s0(:));
        si1 = sum(s1(:));
        H(k) = (u0 - u1)^2/(si0 + si1);     
        k = k + 1;
    end
    kopt = find(H == max(H));
    T = mean(kopt + mi);
end