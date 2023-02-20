function cp=hankel2vec(mtx,han_num)
%% Converting a Hankel into a vector
% Author : Xiaobo Qu
% Email  : quxiaobo@xmu.edu.cn
% Created in 2015
% References: [1] Xiaobo Qu*, Maxim Mayzel, Jian-Feng Cai, Zhong Chen, Vladislav Orekhov*. Accelerated NMR spectroscopy with low-rank reconstruction, Angewandte Chemie International Edition, 54(3):852-854, 2015.
    len_vec=size(mtx,2)+han_num-1;
    cp=zeros(len_vec,1);
    for iter=1:size(mtx,2)
        idx=iter:iter+han_num-1;
        cp(idx)=cp(idx)+mtx(:,iter);        
    end
end
