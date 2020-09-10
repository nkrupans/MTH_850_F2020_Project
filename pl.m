a = dir('u0*.txt');
b = dir('v0*.txt');

u = load(a(1).name);   
n = length(u(:,1))-1;

% Assumes [x,y] \in [0,1]^2
h = 1/n;
[x,y] = meshgrid((0:h:1));

for i = 1:length(a)
    u = load(a(i).name);   
    v = load(b(i).name);
    w = -(v(2:end-1,3:end)-v(2:end-1,1:end-2))/(2*h) ...
        + (u(3:end,2:end-1)-u(1:end-2,2:end-1))/(2*h);
    surf(x(2:end-1,2:end-1),y(2:end-1,2:end-1),w)
    shading interp
    view(0,90)
    %contour(x(2:end-1,2:end-1),y(2:end-1,2:end-1),w,linspace(-10,10,100))
    %colorbar
    %axis equal 
    %mesh(x,y,sqrt(u.^2+v.^2))
    colorbar
    axis equal
    colorbar
    caxis([0 10])
    drawnow
    % uncomment if you want to generate frames for movies 
    %
    % print(['frame' num2str(i,'%05d') '.jpg'],'-djpeg')
    % Once you have generated the jpg files you can compile them
    % into a movie by  
    % ffmpeg -r 10 -i frame%05d.jpg test.mp4
    % 
end
    