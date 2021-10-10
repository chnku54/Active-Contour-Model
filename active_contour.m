function my_active_contour(img)

if(size(img,3) == 3) 
    img=rgb2gray(img); 
end

numPoints = 100;
splinePoints = get_initial_contour(img, numPoints);
X_spline  = splinePoints(1,:);
Y_spline  = splinePoints(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computer the energy term(s) for the snake model

alpha = 0.1; 
beta = 0.005; 
delta = 1; 

points = splinePoints;

x = X_spline'; 
y = Y_spline';

n = length(x);

sigma=1;

 a = [-1 0 1 ; -1 0 1; -1 0 1];
 b = [-1 -1 -1 ; 0 0 0; 1 1 1];

 gau_f = exp(-(a.^2+b.^2)/(2*sigma^2));
 gau_f = gau_f./sum(gau_f,'all');

[p,q] = size(img);

img_gau = zeros(p,q);

row = zeros(1,q)+255;
img_rim = [row; img; row];
column = zeros(p+2,1)+255;
img_rim = [column img_rim column];

img_rim = double(img_rim);

for i = 2:p+1
    for j = 2:q+1
        img_gau(i-1,j-1) = round(sum(img_rim(i-1:i+1,j-1:j+1).*gau_f,'all'));
    end
end

img_gau = double(imfilter(img, gau_f, 'replicate'));

robert_x = [-1 0; 0 1];
robert_y = [0 -1; 1 0];

Gradient_x = zeros(p,q);
Gradient_y = zeros(p,q);

for i = 1:p-1
    for j = 1:q-1
        Gradient_x(i,j) = img(i+1,j+1) - img(i,j);
        Gradient_y(i,j) = img(i+1,j) - img(i,j+1);
    end
end

Gradient_mag = -sqrt(Gradient_x.^2+Gradient_y.^2);

Gradient_mag_max = max(max(Gradient_mag)); 
Gradient_mag_min = min(min(Gradient_mag));

Gradient_mag = (Gradient_mag - Gradient_mag_min)/(Gradient_mag_max - Gradient_mag_min) * (0 - (-1)) + (-1);

Gradient2_x = zeros(p,q);
Gradient2_y = zeros(p,q);

for i = 1:p-1
    for j = 1:q-1
        Gradient2_x(i,j) = Gradient_mag(i+1,j+1) - Gradient_mag(i,j);
        Gradient2_y(i,j) = Gradient_mag(i+1,j) - Gradient_mag(i,j+1);
    end
end

points(3,:) = alpha; 
points(4,:) = beta;
h = 1;


% Set your own number of iterations
numIterations = 1000;
figure
for i=1:numIterations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the locations of the contour
    % You are required to implement the interpolation algorithm by yourself
    % (The nearest interpolation are allowed to used.)
    
        alpha_vec = zeros(1,length(x))+alpha;
        beta_vec = zeros(1,length(x))+beta;
        
        coeff1 = beta_vec/h^2; 
        coeff2 = - alpha_vec/h -4*beta_vec/h^2 ; 
        coeff3 = 2*alpha_vec/h + (6*beta_vec)/h^2; 
        coeff4 = - alpha_vec/h -4*beta_vec/h^2 ; 
        coeff5 = beta_vec/h;

        mat_A = diag(coeff4(1),-(n-1))+diag(coeff5(1:2),-(n-2))+diag(coeff1(1:n-2),-2)+diag(coeff2(1:n-1),-1)+diag(coeff3)+diag(coeff4(1:n-1),1)+diag(coeff5(1:n-2),2)+diag(coeff1(1:2),n-2)+diag(coeff2(1),n-1);

        new_ex = zeros(1,n);
        new_ey = zeros(1,n);

        %Bilinear Interpolation
        for iter = 1:n
          iter_y = x(iter);
          lower_y = floor(iter_y);
          upper_y = ceil(iter_y);
          iter_x = y(iter);
          lower_x = floor(iter_x);
          upper_x = ceil(iter_x);
          new_ex(iter) = (upper_y-iter_y)*((upper_x-iter_x)*Gradient2_x(lower_x,upper_y)+(iter_x-lower_x)*Gradient2_x(upper_x,upper_y))...
              + (iter_y-lower_y)*((upper_x-iter_x)*Gradient2_x(lower_x,lower_y)+(iter_x-lower_x)*Gradient2_x(upper_x,lower_y))
          new_ey(iter) = (upper_y-iter_y)*((upper_x-iter_x)*Gradient2_y(lower_x,upper_y)+(iter_x-lower_x)*Gradient2_y(upper_x,upper_y))...
              + (iter_y-lower_y)*((upper_x-iter_x)*Gradient2_y(lower_x,lower_y)+(iter_x-lower_x)*Gradient2_y(upper_x,lower_y))
        end

        new_ex = new_ex';
        new_ey = new_ey';

        x = inv(mat_A + eye(n)/delta) * (x/delta + new_ex); 
        y = inv(mat_A + eye(n)/delta) * (y/delta + new_ey);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    imshow(img); 
    hold on;
    plot([x; x(1)], [y; y(1)], 'r-');
    hold off;
    pause(0.001)    
end

%{
imshow(I); 
hold on;
plot([xs; xs(1)], [ys; ys(1)], 'r-');
hold off;
%}