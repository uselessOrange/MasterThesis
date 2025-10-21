function locs=compressorCycles(X)
n=length(X);
previous_sample=0;
count=1;
for i = 1:n
    current_sample=X(i);
    if previous_sample == 0 && current_sample ==1
    locs(count)=i;
    count=count+1;
    end
    previous_sample=current_sample;
end
end