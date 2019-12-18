function [] = wrt_txt_orbo(t,v,dv,name,path)



num =length(t);

FID = fopen(fullfile(path,[name '_ORBODATA'  '.orb']),'w+');
fprintf (FID,'STAR: %s \n',name);

for i = 1: num
fprintf (FID,'A %f %f %f\n' ,t(i),v(i),dv(i));
end

fprintf (FID,'END');

fclose(FID);

