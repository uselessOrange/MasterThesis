function x = FeatureExtract(data,locs)

n_locs=length(locs);
f2=zeros(n_locs-1,1);
f1=zeros(n_locs-1,2);
for i=1:n_locs-1
    f1(i,:)=[max(data(locs(i):locs(i+1))),min(data(locs(i):locs(i+1)))];
    f2(i)=mean(data(locs(i):locs(i+1)));
end

x=[f1(:,1)';f1(:,2)';f2'];

end