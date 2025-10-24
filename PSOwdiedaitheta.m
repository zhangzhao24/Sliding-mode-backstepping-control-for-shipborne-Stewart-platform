%% 初始化种群  
clear
%% Sphere
clear
N = 100;                         % 初始种群个数  
d = 2;                          % 可行解维数  
ger = 100;                      % 最大迭代次数       
limit = [0,100];               % 设置alpha位置参数限制  
limit1 = [1e-3,100];               % 设置sigma位置参数限制  
%limit1 = [1e-10,1];               % 设置sigma位置参数限制  
vlimit = [-10, 10];               % 设置速度限制  
w = 0.9;                        % 惯性权重  
wmax=0.9;
wmin=0.4;
c1 = 1;                       % 自我学习因子  
c2 = 1;                       % 群体学习因子   
%x = limit(1) + (  limit( 2 ) -  limit( 1)  ) .* rand(N, d);%初始种群的位置
x(:,1) = limit(1) + (  limit( 2 ) -  limit( 1)  ) .* rand(N, 1);%初始种群的位置
x(:,2) = limit1(1) + (  limit1( 2 ) -  limit1( 1)  ) .* rand(N, 1);%初始种群的位置
v = rand(N, d);                  % 初始种群的速度  
xm = x;                          % 每个个体的历史最佳位置  
ym = zeros(1, d);                % 种群的历史最佳位置  
fxm = ones(N, 1)*inf;               % 每个个体的历史最佳适应度   
fym = inf;                       % 种群历史最佳适应度  
%% 群体更新  
iter = 1;  
% record = zeros(ger, 1);          % 记录器  
figure(1)
while iter <= ger  
  fx = ones(N, 1)*inf;
     %fx = f( x(:,1),x(:,2) ) ;% 个体当前适应度     
     for i = 1:N   
         fprintf("第%d个",i)
              fx(i) = shiyingduJtheta( x(i,1),x(i,2) ) ;% 个体当前适应度   
        if  fx(i)  <fxm(i) 
            fxm(i) = fx(i);     % 更新个体历史最佳适应度  
            xm(i,:) = x(i,:);   % 更新个体历史最佳位置(取值)  
        end   
     end  
    if   min(fxm)<  fym
        [fym, nmin] = min(fxm);   % 更新群体历史最佳适应度  
        ym = xm(nmin, :);      % 更新群体历史最佳位置  
    end  
    w=wmax-(wmax-wmin)*(iter/ger)*(iter/ger);
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% 速度更新 
     % 边界速度处理  
    v(v > vlimit(2)) = vlimit(2);  
    v(v < vlimit(1)) = vlimit(1); 
  
    x = x + v;% 位置更新  
    for i = 1:N
    % 边界位置处理  
    %x(x > limit(2)) = limit(2);  
    %x(x < limit(1)) = limit(1);  
        if (x(i,1)>limit(2))
        x(i,1)= limit(2);
    end
    if (x(i,1)< limit(1))
        x(i,1)= limit(1);
    end
    if (x(i,2)>limit1(2))
        x(i,2)= limit1(2);
    end
    if (x(i,2)< limit1(1))
        x(i,2)= limit1(1);
    end
    end
        record(iter) = fym;%最大值记录 
            plot(record);title('最优适应度变化过程')  
                pause(0.01)  
        iter = iter+1; 
end
    disp(['最优值：',num2str(fym)]);  
%disp(['变量alpha和sigma取值：',num2str(ym)]);
disp(['变量alpha和theta取值：',num2str(ym)]);
%输出结果
alpha = ym(1);
theta = ym(2);
[KL,L]=findLMItheta(alpha,theta);