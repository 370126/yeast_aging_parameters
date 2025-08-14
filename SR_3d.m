% %% 设置导入选项并导入数据
% opts = delimitedTextImportOptions("NumVariables", 4);
% 
% % 指定范围和分隔符
% opts.DataLines = [2, Inf];
% opts.Delimiter = ",";
% 
% % 指定列名称和类型
% opts.VariableNames = ["gene", "mean_rela", "steepness_rela", "skewness_rela"];
% opts.VariableTypes = ["string", "double", "double", "double"];
% 
% % 指定文件级属性
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % 指定变量属性
% opts = setvaropts(opts, "gene", "WhitespaceRule", "preserve");
% opts = setvaropts(opts, "gene", "EmptyFieldRule", "auto");
% 
% % 导入数据
% lifes_set_alpha_bp_short_df = readtable("C:\Users\20145\Desktop\westlake_summer\lifes_set_alpha_bp_short_df.csv", opts);
% 
% opts = delimitedTextImportOptions("NumVariables", 5);
% 
% % 指定范围和分隔符
% opts.DataLines = [2, Inf];
% opts.Delimiter = ",";
% 
% % 指定列名称和类型
% opts.VariableNames = ["para", "eff", "lifespan_rela", "steepness_rela", "skewness_rela"];
% opts.VariableTypes = ["categorical", "categorical", "double", "double", "double"];
% 
% % 指定文件级属性
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % 指定变量属性
% opts = setvaropts(opts, ["para", "eff"], "EmptyFieldRule", "auto");
% 
% % 导入数据
% simu_long = readtable("C:\Users\20145\Desktop\westlake_summer\simu_long.csv", opts);
% 
% clear opts
% 

%%
load('mat_data.mat')

%% plot

s=5*ones(height(lifes_set_alpha_bp_short_df),1);
figure
for i=1:length(experi_para)
    y_plot=steepness_rela(i,:);
    z_plot=skewness_rela(i,:);
    plot(y_plot,z_plot,'LineWidth',2)
    hold on
end
title('Relative Skewness v.s. Relative Steepness')
xlabel('steepness_rela')
ylabel('skewness_rela')
legend(experi_para)
scatter(lifes_set_alpha_bp_short_df.steepness_rela, ...
    lifes_set_alpha_bp_short_df.skewness_rela, ...
    s,"filled", ...
    'MarkerFaceAlpha',0.5)
xlabel('steepness')
ylabel('skewness')
axis on

% figure
% for i=1:length(experi_para)
%     x_plot=lifespan_rela(i,:);
%     y_plot=steepness_rela(i,:);
%     z_plot=skewness_rela(i,:);
%     plot3(x_plot,y_plot,z_plot,'LineWidth',2)
%     hold on
% end
% title('Relative Steepness, Relative Lifespan, Relative Skewness')
% xlabel('lifespan_rela')
% ylabel('steepness_rela')
% zlabel('skewness_rela')
% legend(experi_para)
% scatter3(lifes_set_alpha_bp_short_df.mean_rela, ...
%     lifes_set_alpha_bp_short_df.steepness_rela, ...
%     lifes_set_alpha_bp_short_df.skewness_rela, ...
%     s,"filled", ...
%     'MarkerFaceAlpha',0.5)
% xlabel('lifespan')
% ylabel('steepness')
% zlabel('skewness')
% axis on
