function side=side_of_flat(p1,p0,ori)
side=sign(sum(repmat(ori,size(p1,1),1).*(p1-p0),2));
end