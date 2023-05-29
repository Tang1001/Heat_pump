function [aver_date,aver_val]=aver_min2hour(date,val)

index_aver=1;
aver_date(index_aver)=date(1);
aver=val(1);
num_aver=1;
for i=1:length(date)-1
    if hour(date(i))~=hour(date(i+1))
        aver_val(index_aver)=aver/num_aver;
        aver_date(index_aver+1)=date(i+1);
        index_aver=index_aver+1;
        aver=val(i+1);
        num_aver=1;
    else
        aver=aver+val(i+1);
        num_aver=num_aver+1;
    end

    %The end
    if i==length(date)-1
        aver_val(index_aver)=aver/num_aver;
    end

end

