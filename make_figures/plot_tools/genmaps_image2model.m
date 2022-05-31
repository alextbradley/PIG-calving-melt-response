function [model2image, image2model] = genmaps_image2model(ximage, yimage, xmod, ymod)
Rimage = [ximage'; yimage'];
Rmodel = [xmod'; ymod'];
model_length = sqrt(diff(Rmodel(1,:)).^2 + diff(Rmodel(2,:))^2);
image_length = sqrt(diff(Rimage(1,:)).^2 + diff(Rimage(2,:))^2);
%image2model

shift = Rimage(:,1);
model_angle = atand(diff(ymod)/diff(xmod));
image_angle = atand(diff(yimage)/diff(ximage));
angle = model_angle - image_angle + 180;
rotation = [cosd(angle), -sind(angle);   
            sind(angle), cosd(angle)];
        

reflection = [cosd(2*model_angle), sind(2*model_angle);
            sind(2*model_angle), -cosd(2*model_angle)];

       
scaling = model_length/image_length;

second_shift = Rmodel(:,1);

image2model = @(u)  scaling*reflection*rotation*(u - shift) + second_shift;
model2image = @(U) scaling^(-1) * inv(rotation) * inv(reflection) * (U - second_shift) +shift;
end
