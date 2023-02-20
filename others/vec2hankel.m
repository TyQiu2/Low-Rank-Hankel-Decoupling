function cp=vec2hankel(vec,han_num)
%% Converting a vector into a Hankel matrix
% Author : Xiaobo Qu
% Email  : quxiaobo@xmu.edu.cn
% Created in 2015
% References: [1] Xiaobo Qu*, Maxim Mayzel, Jian-Feng Cai, Zhong Chen, Vladislav Orekhov*. Accelerated NMR spectroscopy with low-rank reconstruction, Angewandte Chemie International Edition, 54(3):852-854, 2015.
%%
if size(vec,1)==1
    vec=reshape(vec,[],1);
end
cp=[];
len_vec=length(vec);
    for iter=1:len_vec-han_num+1
        vec_one_row=vec(iter:iter+han_num-1);
        cp=[cp vec_one_row];
    end
end
%%
