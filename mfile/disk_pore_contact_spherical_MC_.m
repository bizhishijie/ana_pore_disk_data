rc_p=zeros(3,num_point_all);
idx_contact=zeros(2,num_point_all);
rs=zeros(1,num_point_all);

range=[min(Rc_in,[],2) max(Rc_in,[],2)];

tic
parfor ii=1:num_point_all
    point_in_pore=0;
    while point_in_pore==0
        rc_p_test=rand(3,1).*[range(:,2)-range(:,1)]+range(:,1);
        if inside_polyhedron_(rc_p_test,Rc_in)
            rc2=rc_p_test-Rc;
            d_z=dot(rc2,Ori);
            rc_xy=rc2-repmat(d_z,3,1).*Ori;
            d_xy=sqrt(sum(rc_xy.^2,1));

            if any(abs(d_z)<H/2&d_xy<D/2)
                continue
            else
                idx_face=abs(d_z)>H/2&d_xy<=D/2;
                idx_edge=abs(d_z)>H/2&d_xy>D/2;
                idx_sidewall=abs(d_z)<=H/2&d_xy>D/2;

                rs_all=zeros(1,size(Rc,2));
                rs_all(idx_face)=abs(d_z(idx_face))-H/2;
                rs_all(idx_edge)=sqrt((d_z(idx_edge)-H/2).^2+(d_xy(idx_edge)-D/2).^2);
                rs_all(idx_sidewall)=d_xy(idx_sidewall)-D/2;

                point_in_pore=1;
            end
        else
            continue
        end
    end

    rc_p(:,ii)=rc_p_test;
    [rs_tmp,idx]=min(rs_all);
    rs(ii)=rs_tmp;
    idx_contact(:,ii)=[idx_face(idx)*1+idx_edge(idx)*2+idx_sidewall(idx)*3;idx];
end
toc

