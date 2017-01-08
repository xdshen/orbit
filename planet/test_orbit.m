function  test_orbit2

%% һ���趨���ǲ���
% planet: image,r,R,T_revolve,T_spin, sunID,

% Sun
data{1}=struct( 'name','sun', ...
                'image', 'sun.jpg',...
                'radius', 6.955e5,...
                'L',0,...
                'T_revolve',inf,...
                'T_spin',25,...
                'sunID',0,...
                'mass',1.9891e30);
% mercury
data{2}=struct( 'name','mercury', ...
                'image', 'mercury.jpg',...
                'radius', 2440,...
                'L',5.7909100e7 ,...      % km
                'T_revolve',87.969,...
                'T_spin',58.646,...
                'sunID',1,...
                'mass',3.3022e23);            
            
%  venus
data{3}=struct('name','venus', ...
                'image', 'venus.jpg',...
                'radius', 2440,...
                'L',1.08208000e8 ,...      % km
                'T_revolve',243,...
                'T_spin',-224.7,...
                'sunID',1,...
                'mass',4.8675e24);           
% earth
data{4}=struct('name','earth', ...
                'image', 'earth.jpg',...
                'radius', 6378,...
                'L',1.4960e+08 ,...      % km
                'T_revolve',365.24219,...
                'T_spin',1,...
                'sunID',1,...
                'mass',5.965e24);           
% moon
data{5}=struct('name','moon', ...
                'image', 'moon.jpg',...
                'radius', 1738,...
                'L',384400 ,...      % km
                'T_revolve',27.32,...
                'T_spin',27.32,...
                'sunID',4,...
                'mass',7.349e22);    
% mars
data{6}=struct('name','mars', ...
                'image', 'mars.jpg',...
                'radius', 6794/2,...
                'L',2.2793664e8 ,...      % km
                'T_revolve',687,...
                'T_spin',24.6229/24,...
                'sunID',1,...
                'mass',6.4219e23);    
% % jupiter
% data{7}=struct('name','jupiter', ...
%                 'image', 'jupiter.jpg',...
%                 'radius', 142987/2,...
%                 'L',5.20336301*1.4960e+08 ,...      % km
%                 'T_revolve',11.86*365,...
%                 'T_spin',(9+50/60)/24,...
%                 'sunID',1,...
%                 'mass',1.9e27);    
% % saturn
% data{8}=struct('name','saturn', ...
%                 'image', 'saturn.jpg',...
%                 'radius', 120540/2,...
%                 'L',1.433449369500000e+09 ,...      % km
%                 'T_revolve',29.53216*365,...
%                 'T_spin',10.546/24,...
%                 'sunID',1,...
%                 'mass',5.688e26);    

% ################## Ϊ����ʾ�ÿ�

    for i=1:length(data)
        switch data{i}.name
            case 'sun'
                % ̫���İ뾶��С70��
                data{i}.radius = data{i}.radius/70;
            case {'mercury','venus','earth','mars'}
                % �����������⣬���ھ�����С4000��
                data{i}.L = data{i}.L/4000;
            case 'moon'
                % ���µľ�����С40��
                data{i}.L = data{i}.L/40;
            case 'jupiter'
                % Jupiter̫����
                data{i}.L = data{i}.L/4000;
                data{i}.radius = data{i}.radius/10;                
            case 'saturn'
                % Saturn̫����
                data{i}.L = data{i}.L/4000;
                data{i}.radius = data{i}.radius/10;        
                
            otherwise
        end
    end
    

% ###############################          
            
            
%% ������ʼ��   
N=50;
n_step=0;


for i=1:length(data)
   % load image
   data{i}.I = load_image(data{i}.image);
   % create sphere
   data{i}.sphere = create_sphere(data{i}.radius,N); 
end          

% ��ʼ��figure
h = figure();
set(h,'outerposition',get(0,'screensize'));
grid on;
set(h,'Resize','off');
L_max = data{length(data)}.L;      % ��Զ����ʾ�ľ���(����һ���Ǿ�����Զ�ģ�)
axis([-L_max L_max -L_max L_max -L_max L_max]);    
view([L_max -L_max L_max]);
hold on;
% ��ʼ��surf
for i=1:length(data)
    data{i}.hSurf = surf(data{i}.sphere.x,data{i}.sphere.y,data{i}.sphere.z,data{i}.I,'FaceColor','texturemap','edgecolor','none');    
end



%% ������ʼѭ����ͼ
for t = 1:0.1:30000   
    n_step=n_step+1;
   
    % 1��ѡȷ��ÿ�����ǵ�����λ��
    c = zeros(length(data),2);
    for i=1:length(data)
      if data{i}.sunID==0
          c(i,:)=[0 0];
      else
          j = data{i}.sunID;
          c(i,:) = c(j,:) + data{i}.L* [cos(t/data{i}.T_revolve), sin(t/data{i}.T_revolve)];
      end
      data{i}.pos = c(i,:);
    end
    
    % 2��ȷ������mesh/surf�����飻
    p = cell(length(data),1);    
    for i=1:length(data)
       
       data{i}.sphere2 =  transform(data{i}.sphere, c(i,:), t/data{i}.T_spin);          % �µ�sphere
    end
    
    % 3����̬��ͼ     

    data = drawimage_texture(data);
    data = plot_trace(data);
 
    
    %axis equal;
    drawnow;

end
hold off;  
  
end
%%


function data = drawimage_texture(data)
% ��ͼƬ������
    for i = 1:length(data)
    	if ~isfield(data{i},'hSurf')
            data{i}.hSurf = surf(data{i}.sphere2.x,data{i}.sphere2.y,data{i}.sphere2.z,data{i}.I,'FaceColor','texturemap','edgecolor','none');    
        else
            set(data{i}.hSurf,'xdata',data{i}.sphere2.x,'ydata',data{i}.sphere2.y,'zdata',data{i}.sphere2.z,'cdata',data{i}.I,'FaceColor','texturemap','edgecolor','none');
        end
    end
    
end

function I = load_image(filepath)
% ����ͼƬ
    I = imread(['images\' filepath]);
    I = flipdim(I, 1);
end
  function s = create_sphere(R, N)
  % ����һ����
        [x,y,z]=sphere(N);
        s.x = R*x;
        s.y = R*y;
        s.z = R*z;
  end
  
  function p = transform(s, central, theta)
  % ��������λ�ú���ת�ĽǶȷ��ر����ļ��ϣ�����һ��mesh/surf
    t = rotate_z(s, theta);
    p.x = t.x + central(1);
    p.y = t.y + central(2);
    p.z = t.z;
  end
  
  function p = rotate_z(s, theta)
  % theta in radians, ����Z��תtheta�Ƕ�
    t = (s.x + s.y * 1i) * exp(1i*theta);
    p.x = real(t);
    p.y = imag(t);
    p.z = s.z;
  end
  
  function data = plot_trace(data)
  % �����е����ǵĹ켣
    for i =1:length(data)
        if ~(data{i}.sunID==0)
            plot_trace_i(i,data{i}.L,data{data{i}.sunID}.pos,0);                
        end            
    end
      function plot_trace_i(i,L,c_xyz,theta)
          u=(0:0.1:2*pi).'+theta;
          [c_xy2, u2]= meshgrid(c_xyz(1:2),u);
          T = repmat([L,L],length(u),1);
          pos = c_xy2+T.*[cos(u2(:,1)) sin(u2(:,2))];
          if ~isfield(data{i},'hPlot3')
              data{i}.hPlot3 = plot3(pos(:,1),pos(:,2),zeros(size(u)),'--');
          else
              set(data{i}.hPlot3,'xdata',pos(:,1),'ydata',pos(:,2),'zdata',zeros(size(u)));
          end
      end
  end