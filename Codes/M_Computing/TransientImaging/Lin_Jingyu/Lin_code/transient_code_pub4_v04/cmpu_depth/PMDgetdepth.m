% PDM depth estimation

sel_data = '1';
iter = 20; % measuring times

w = 165;
h = 120;
count = w*h;
pmddepth = zeros(h,w-2,iter);
for i=1:iter
    % acquiring depth
%     cmdpath = '..\\getDepth\\Release';
    cmdpath = '..\TOF_tools';
    captool = 'getDepth.exe';
    syscmd = sprintf('%s\\%s',cmdpath,captool);
    system(syscmd);

    % reading data
    fn = 'depth.dat';
    [fid, msg] = fopen(fn,'r');
    if fid==-1
        fprintf(msg);
        fprintf('\n');
        dat = 0;
        return;
    end

    fdat = fread(fid, count, 'float32');
    fclose(fid);
    depthmap = reshape(fdat,[w h]);
    pmddepth(:,:,i) = depthmap(1:163,end:-1:1).';

    % show
    figure(1); imagesc(pmddepth(:,:,i)); colorbar
%     figure(1); imagesc(pmddepth(40:80,80:100,i)); colorbar
    title(i)
    drawnow
end

fn = sprintf('results\\pmddepth_%s',sel_data);
save(fn,'pmddepth');
return
