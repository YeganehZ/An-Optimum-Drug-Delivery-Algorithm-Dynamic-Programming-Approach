
clear all
D=0.43;           % Diffussion rate
M=10^6;           % Derug rate 
T=1; 
N=25;           % Number of voxels
do=2;           % distante to organ 
dt=1;           % distance to tumor
Q=20000;      % thershould
% Calaulate yo and yt for the input as a one step
yo=zeros(1,N);
yt=zeros(1,N);
yo(1)=(M/(4*pi*D*do))*erf(sqrt(do^2/(4*D*T)));
yt(1)=(M/(4*pi*D*dt))*erf(sqrt(dt^2/(4*D*T)));
for m=2:N
    yo(m)=(M/(4*pi*D*do))*(erf(sqrt(do^2/(4*D*(m-1)*T)))-erf(sqrt(do^2/(4*D*m*T))));
    yt(m)=(M/(4*pi*D*dt))*(erf(sqrt(dt^2/(4*D*(m-1)*T)))-erf(sqrt(dt^2/(4*D*m*T))));
end
%yo
%yt

% Intialization: k=1
WO1(1)=yo(1);
WT1(1)=yt(1);
WO2(1)=yo(1);
WT2(1)=yt(1);
X1(1)=1;
X2(1)=1;

% After k>1

for k=2:N
    
    X11=[X1,1];
    X12=[X1,0];
    X21=[X2,1];
    X22=[X2,0];
    
    BO11=0;
    BO12=0;
    BO21=0;
    BO22=0;
    BT11=0;
    BT12=0;
    BT21=0;
    BT22=0;
    
    % Calculate branch metrics
    for m=1:k
        
        BO11=BO11+X11(m)*yo(k+1-m);
        BO12=BO12+X12(m)*yo(k+1-m);
        BO21=BO21+X21(m)*yo(k+1-m);
        BO22=BO22+X22(m)*yo(k+1-m);
        
        BT11=BT11+X11(m)*yt(k+1-m);
        BT12=BT12+X12(m)*yt(k+1-m);
        BT21=BT21+X21(m)*yt(k+1-m);
        BT22=BT22+X22(m)*yt(k+1-m);
        
    end
    
    % Calculate transition weight
    W11=WO1(k-1)+BO11;
    W12=WO1(k-1)+BO12;
    W21=WO2(k-1)+BO21;
    W22=WO2(k-1)+BO22;
    
    % Select survival paths
    
    % Survival path for "Injection" stste or state "1"
    
    if  W21<W11
        
        if  BT21>=Q
            X1=X21;
            WO1(k)=W21;
            WT1(k)=BT21;
        
        elseif BT11>=Q           
            X1=X11;
            WO1(k)=W11;
            WT1(k)=BT11;    
            
        elseif BT11>BT21           
            X1=X11;
            WO1(k)=W11;
            WT1(k)=BT11;
            
        else
            X1=X21;
            WO1(k)=W21;
            WT1(k)=BT21;
        end
        
    else
        if  BT11>=Q
            X1=X11;
            WO1(k)=W11;
            WT1(k)=BT11;
            
        elseif BT21>=Q
            X1=X21;
            WO1(k)=W21;
            WT1(k)=BT21;
            
        elseif BT11>BT21           
            X1=X11;
            WO1(k)=W11;
            WT1(k)=BT11;
            
        else
            X1=X21;
            WO1(k)=W21;
            WT1(k)=BT21;
        end
    end
    
    % Survival path for "Stop-Injection" stste or state "2"
    
     if  W22<W12
        
        if  BT22>=Q
            X2=X22;
            WO2(k)=W22;
            WT2(k)=BT22;
        
        elseif BT12>=Q           
            X2=X12;
            WO2(k)=W12;
            WT2(k)=BT12;    
            
        elseif BT12>BT22           
            X2=X12;
            WO2(k)=W12;
            WT2(k)=BT12;
            
        else
            X2=X22;
            WO2(k)=W22;
            WT2(k)=BT22;
        end  
        
    else
        if  BT12>=Q
            X2=X12;
            WO2(k)=W12;
            WT2(k)=BT12;
            
        elseif BT22>=Q
            X2=X22;
            WO2(k)=W22;
            WT2(k)=BT22;
            
        elseif BT12>BT22           
            X2=X12;
            WO2(k)=W12;
            WT2(k)=BT12;
            
        else
            X2=X22;
            WO2(k)=W22;
            WT2(k)=BT22;
        end
    end                
end

% Final decision

if  WO1(N)<WO2(N)
    if  WT1(N)>=Q
        X=X1;
        WO=WO1;
        WT=WT1;
        
    elseif WT2(N)>=Q
        X=X2;
        WO=WO2;
        WT=WT2;
        
    elseif  WT1(N)>WT2(N)
        X=X1;
        WO=WO1;
        WT=WT1;
    else
        X=X2;
        WO=WO2;
        WT=WT2;
    end
  
else
    if  WT2(N)>=Q
        X=X2;
        WO=WO2;
        WT=WT2;
        
    elseif WT1(N)>=Q
        X=X1;
        WO=WO1;
        WT=WT1;
        
    elseif  WT1(N)>WT2(N)
        X=X1;
        WO=WO1;
        WT=WT1;
    else
        X=X2;
        WO=WO2;
        WT=WT2;
    end
end

n=1:N;

plot(n,WO,'b')
xlabel('Delivery Duration (NT)'); ylabel('Sum of Absorbed Molecules by Organ');
hold on

figure
plot(n,WT,'b')
xlabel('Delivery Duration (NT)'); ylabel('Number of Molecules at Tumor');
hold on
figure  
for m=1:N
    for i=1:10
        XX((m-1)*10+i)=M*X(m);
    end
end
m=1:10*N;
plot(m,XX,'b')
xlabel('Delivery Duration (NT)'); ylabel('Level of Injection Molecules ');        
            
            
            
        
        
    
    
    
    
    
    
    
    
        
        
        
        
        
        
        
        
        
        
      

 
 