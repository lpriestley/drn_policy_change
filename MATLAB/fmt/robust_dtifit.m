function robust_dtifit(dirname)
%
% robust_dtifit(dirname)
%
if exist('~/fsl/src/etc') == 7,
   addpath('~/fsl/src/etc/matlab');
else
   addpath('/usr/local/fsl/etc/matlab')
end

datafile = [dirname '/data'];
maskfile = [dirname '/nodif_brain_mask'];
bvecsfile = [dirname '/bvecs'];
bvalsfile = [dirname '/bvals'];


% Check format of data and read mask and some other stuff
tmpfname = prepare_read_avw_img_slice(datafile);
[dims,scales] = read_avw_hdr(tmpfname);
mask = read_avw(maskfile);
if size(mask,1)~=dims(1) || size(mask,2)~=dims(2) || size(mask,3)~=dims(3),
   error('Size mismatch between mask and data');
end
bvecs = load(bvecsfile,'-ascii');
bvals = load(bvalsfile,'-ascii');
bvecs = bvecs';
bvals = bvals';

% allocate
s0 = zeros(size(mask));
l1 = zeros(size(mask));
l2 = zeros(size(mask));
l3 = zeros(size(mask));
v1 = zeros([size(mask) 3]);
v2 = zeros([size(mask) 3]);
v3 = zeros([size(mask) 3]);


disp('running robust dtifit...');
% run
b0indx = find(bvals==0);
dwindx = find(bvals~=0);

pr = [-7.29 4; -7.30 4; -7.31 4];
pri = [6 7 8];


for sl=1:dims(3)
  maskind = find(mask(:,:,sl));
  [xind,yind] = ind2sub(size(mask(:,:,sl)),maskind);
  sldata = read_avw_img_slice(tmpfname,sl);
  for i = 1:length(maskind)    
    y = abs(squeeze(sldata(xind(i),yind(i),:)));
    op = [log(var(y(dwindx))); [log(abs(mean(y(b0indx)))) 0 0 0 log([0.00075 0.00070 0.00065])]'];
    tp = df_MAP([0.2*op(1); op(2:end)],y,bvals,bvecs,1,[0.2*op(1) 0.00001; op(2) 0.00001; pr],[1 2 pri]);
    op = [op(1); tp(2:end)];
    
    [np,ney] = df_MAP(op,y,bvals,bvecs,1,pr,pri);
    [rp,rey] = df_MAP(np,y,bvals,bvecs,1,pr,pri);
    
    
    s0(xind(i),yind(i),sl) = exp(rp(2));

        
    % sort eigenvalues
    [rl,il] = sort(exp(rp(6:8)));
    
    
    R = phi2R(rp(3:5));
    v1(xind(i),yind(i),sl,:) = R(:,il(3))'; %(R*[1;0;0])';
    v2(xind(i),yind(i),sl,:) = R(:,il(2))'; %(R*[0;1;0])';
    v3(xind(i),yind(i),sl,:) = R(:,il(1))'; %(R*[0;0;1])';
    
    l1(xind(i),yind(i),sl) = rl(3);
    l2(xind(i),yind(i),sl) = rl(2);
    l3(xind(i),yind(i),sl) = rl(1);
  end  
end

% Remove temporary files
cleanup_read_avw_img_slice(tmpfname);

disp('saving results...');
% saving
save_avw(s0,[dirname '/rrdti_S0'],'f',scales);
save_avw(v1,[dirname '/rrdti_V1'],'f',scales);
save_avw(v2,[dirname '/rrdti_V2'],'f',scales);
save_avw(v3,[dirname '/rrdti_V3'],'f',scales);
save_avw(l1,[dirname '/rrdti_L1'],'f',scales);
save_avw(l2,[dirname '/rrdti_L2'],'f',scales);
save_avw(l3,[dirname '/rrdti_L3'],'f',scales);

return

function R = phi2R(phi)

 s1 = sin(phi(1));
 c1 = cos(phi(1));
 s2 = sin(phi(2));
 c2 = cos(phi(2));
 s3 = sin(phi(3));
 c3 = cos(phi(3));
   
R = [c3*c2 -s3*c1-c3*s2*s1 s3*s1-c3*s2*c1;
     s3*c2 c3*c1-s3*s2*s1 -c3*s1-s3*s2*c1;
     s2 c2*s1 c2*c1];

return

%R = [cos(phi(3)) -sin(phi(3)) 0;sin(phi(3)) cos(phi(3)) 0;0 0 1]...
%    *[cos(phi(2)) 0 -sin(phi(2));0 1 0;sin(phi(2)) 0 cos(phi(2))]...
%    *[1 0 0;0 cos(phi(1)) -sin(phi(1));0 sin(phi(1)) cos(phi(1))];






