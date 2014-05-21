%Making video file .avi, .mp4
%
%Example:
% writerObj = VideoWriter('peaks.avi');
% open(writerObj);
% 
% Z = peaks; surf(Z); 
% axis tight
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
% 
% for k = 1:20 
%    surf(sin(2*pi*k/20)*Z,Z)
%    frame = getframe;
%    writeVideo(writerObj,frame);
% end
% 
% close(writerObj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arg in:
%Str= 'name.mp4' file name;
%tf = final time;
%frame = variable it will be show in the movie
function createMovie(str,tf,frame)

writerObj = VideoWriter(str);
open(writerObj);

axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

for k = 1:tf
   writeVideo(writerObj,frame(k));
end

close(writerObj);

