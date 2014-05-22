%Making video file .avi
%
%Example:
% writerObj = VideoWriter('C:\Users\Convidado.Benedito-HP\Desktop\peaks.avi');
% open(writerObj);
% 
% Z = peaks; surf(Z); 
% axis tight
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
% 
% for k = 1:200 
%    surf(sin(2*pi*k/20)*Z,Z)
%    frame = getframe;
%    writeVideo(writerObj,frame);
% end
% 
% close(writerObj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arg in:
%Str= 'C:\Users\John\Desktop\name' file name;
%tf = final time;
%frame = variable it will be show in the movie
function Draft_createMovie(str,tf,frame)
%cd file
writerObj = VideoWriter(str);
open(writerObj);

axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

for k = 1:tf
   writeVideo(writerObj,frame(k));
end

close(writerObj);

