function s=length_spring(n, l)
    s=[];
    s_eq=l/n; %uniformément distribué
    i=1;
    r=mod(n, 2);
    if r==0
        while i<n
            corr=(0.2*rand(1, 1))*s_eq; %facteur correcteur
            s=[s s_eq-corr];
            s=[s s_eq+corr];
            i=i+2;
        end
    else
        while i<n
            corr=(0.2*rand(1, 1))*s_eq; %facteur correcteur
            s=[s s_eq-corr];
            s=[s s_eq+corr];
            i=i+2;
        end
        s(n)=s_eq;
end

        