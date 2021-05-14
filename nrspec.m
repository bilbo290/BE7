clear 
syms beta
y = 1.4;    
delta = 6*pi/180;
M1m = linspace(1,2,51)
m = length(M1m)
data = zeros(m,9)
for k = 1:m
    M1 = M1m(k);
    func = 2*cot(beta)*((M1^2)*sin(beta)^2 - 1)/(2+ (M1^2)*(y + cos(2*beta))) - tan(delta) ;
    n = 40; % iteration
    x = zeros(n,1); 
    x(1) = 30*pi/180;
    for i = 1:n-1
        x(i+1) = x(i) - subs(func,x(i))/subs(diff(func),x(i));
    end
    beta1 = x(n);
    disp(beta1*180/pi);
    if beta1*180/pi > 90 || beta1*180/pi < 0 || beta1*180/pi > 67.5
         p2top1 = ((((y+1)*((M1^2))/((((y-1)*(M1^2))) +2)))^(y/(y-1)))  * ((y+1)/(2*y*(M1^2) - (y-1)))^(1/(y-1));
        beta1 = 90*pi/180;
    else
        p2top1 = ((((y+1)*((M1^2)*sin(beta1)^2))/((((y-1)*(M1^2)*(sin(beta1)^2)) +2)))^(y/(y-1)))  * ((y+1)/(2*y*(M1^2)*(sin(beta1)^2) - (y-1)))^(1/(y-1));
    end
    r2tor1 = ((y+1)*(M1^2)*(sin(beta1)^2))/(2 +((y-1)*(M1^2)*(sin(beta1)^2)));
    T2toT1 = p2top1/r2tor1;



    M2 = (1/sin(beta1-delta))*sqrt((((M1^2)*sin(beta1)^2 + (2/(y-1))))/((2*y*(M1^2)*(sin(beta1)^2)/(y-1))-1));

    x(1) = 30*pi/180;
    func = 2*cot(beta)*((M2^2)*sin(beta)^2 - 1)/(2+ (M2^2)*(y + cos(2*beta))) - tan(delta) ;
    for i = 1:n-1
        x(i+1) = x(i) - subs(func,x(i))/subs(diff(func),x(i));
    end
    beta2 = x(n);
    disp(beta2*180/pi);
    if abs(beta2*180/pi) > 90 || beta2*180/pi < 0
        p3top2 = ((((y+1)*((M2^2))/((((y-1)*(M2^2))) +2)))^(y/(y-1)))  * ((y+1)/(2*y*(M2^2) - (y-1)))^(1/(y-1));
        beta2 = 90*pi/180;
    else
        p3top2 = ((((y+1)*((M2^2)*sin(beta2)^2))/((((y-1)*(M2^2)*(sin(beta2)^2)) +2)))^(y/(y-1)))  * ((y+1)/(2*y*(M2^2)*(sin(beta2)^2) - (y-1)))^(1/(y-1));
    end
    r3tor2 = ((y+1)*(M2^2)*(sin(beta2)^2))/(2 +((y-1)*(M2^2)*(sin(beta2)^2)));
    T3toT2 = p3top2/r3tor2;
    M3 = (1/sin(beta2-delta))*sqrt((((M2^2)*sin(beta2)^2 + (2/(y-1))))/((2*y*(M2^2)*(sin(beta2)^2)/(y-1))-1));
    p4top3 = ((((y+1)*((M3^2))/((((y-1)*(M3^2))) +2)))^(y/(y-1)))  * ((y+1)/(2*y*(M3^2) - (y-1)))^(1/(y-1));
    eta = p5top4*p4top3*p3top2*p2top1;
    data(k,:) = [M1 beta1*180/pi p2top1 M2 beta2*180/pi p3top2 M3 p4top3 eta];
end
plot(data(:,1),data(:,9))
hold on
nr = 1 - 0.075.*(M1m-1).^1.35
plot(data(:,1),nr)
hold off