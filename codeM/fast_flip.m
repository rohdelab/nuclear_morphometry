function img_array_out = fast_flip(img_array_in)

[ny,nx,nz] = size(img_array_in);

for a=1:nz
    xx=img_array_in(:,:,a);
    d_lr=sum(xx); d_lr=d_lr/sum(d_lr);
    d_ud=sum(xx'); d_ud=d_ud/sum(d_ud);
    m_lr=sum(d_lr.*[1:length(d_lr)])-find(max(d_lr)==d_lr);
    m_ud=sum(d_ud.*[1:length(d_ud)])-find(max(d_ud)==d_ud);
    if m_lr<0
        xx=fliplr(xx);
    end
    if m_ud<0
        xx=flipud(xx);
    end
    img_array_out(:,:,a)=xx;
end


end