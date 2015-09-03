function imshowBW(Im,c,lab)
% show image in Black and White. 
% c: range in grayscale
% lab: labels for axis in cell-array
% 
% By Rasmus Dalgas Kongskov, 20/02/2014, DTU

% image-size
[N1,N2] = size(Im); 

% real or not
re = norm(imag(Im))<(N1*eps);

if nargin<2 || isempty(c)
    if re
        Im = real(Im);
        c = [min(Im(:)),max(Im(:))]; % colorlimit
    else
        cr = [min(real(Im(:))),max(real(Im(:)))]; % colorlimit
        ci = [min(imag(Im(:))),max(imag(Im(:)))]; % colorlimit
    end
end

if N2==1
    N1 = sqrt(N1);
    N2 = N1;
end

if re % real image only
    Im = real(Im);
    imagesc(reshape(Im,N1,N2)); colormap('gray');
    axis image off; caxis(c); colorbar
    
else % imag and real image
    subplot(1,2,1);
    imagesc(reshape(real(Im),N1,N2)); colormap('gray');
    axis image off; caxis(cr); title('Real'); colorbar
    subplot(1,2,2);
    imagesc(reshape(imag(Im),N1,N2)); colormap('gray');
    axis image off; caxis(ci); title('Imag'); colorbar
end

if nargin==3 % including labels
    if re % real image only
        
        Im = real(Im);
        imagesc(reshape(Im,N1,N2)); colormap('gray'); axis image;
        xlabel(lab(1)); ylabel(lab(2)); set(gca,'xtick',[],'ytick',[]);
        caxis(c); colorbar
    else % imag and real image
        
        subplot(1,2,1);
        imagesc(reshape(real(Im),N1,N2)); colormap('gray'); axis image;
        xlabel(lab(1)); ylabel(lab(2)); set(gca,'xtick',[],'ytick',[]);
        caxis(cr); title('Real'); colorbar
        subplot(1,2,2);
        imagesc(reshape(imag(Im),N1,N2)); colormap('gray'); axis image
        xlabel(lab(1)); ylabel(lab(2)); set(gca,'xtick',[],'ytick',[]);
        caxis(ci); title('Imag'); colorbar
    end
end

end