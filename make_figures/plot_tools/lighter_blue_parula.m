function cmap = lighter_blue_parula(N,f)
    %parula colormap shifted towards the lighter end of the spectrum. N is the size of the original matrix, f is the fraction of entries removed
    cmap = parula(N);
    idx = round(f*N);
    cmap = cmap(idx:end,:);
end
