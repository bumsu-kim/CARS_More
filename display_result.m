function dispres = display_result(Results, ~, fparam, ProbType, n)

%% PARAMETERS
%  ------------- INPUT ------------
% Results = cell of structs
% Results{1}, ..., Results{N}

% Each result in Results has the following fields
% 1. sol            solution sequence
% 2. objval_seq     objective function value sequence
% 3. gamma_seq      gamma sequence
% 4. num_iter       number of iterations (total, performed)
% 5. num_queries    number of queries at each iteration (cumulative)
% 6. name           name of the method (used for displaying the result)

% Some results may have the following fields:
% 1. Hseq           Hessian estimate sequence (QN)
% 2. gseq           gradient estimate sequence (QN, FDSA)
% 3. 

%%
displog = true;

dispres = struct;
N = length(Results);
% n = length(Results{1}.sol{1});

fmin = fparam.fmin;

% err = cell(N,1);
legend_str = cell(N,1);
for i=1:N
%     err{i} = zeros(Results{i}.num_iter + 1, 1);
%     for j=1:(Results{i}.num_iter + 1)
%         err{i}(j) = norm(Results{i}.sol{j} - true_sol);
%     end
%     disp(['Number of queries / dim (' , Results{i}.name,') = ', num2str(Results{i}.num_queries(end) / n),' with ',num2str(Results{i}.num_iter) , ' iterations']);
    legend_str{i} = Results{i}.name;
end

% Figures Saving Option
% ProbType = 'CVX';
% ProbType = 'RLR';
if strcmp(ProbType,'CVX')
    filepath = ['E:/Google Drive_nnaavverr/STUDY/Research in Optimization/', ...
                '1. 0th order optimization/Figures_Sept/convex'];
elseif strcmp(ProbType,'RLR')
    filepath = ['E:/Google Drive_nnaavverr/STUDY/Research in Optimization/' ...
                '1. 0th order optimization/Figures_Sept/nonconvex'];
elseif  strcmp(ProbType,'QRT')
    filepath = ['E:/Google Drive_nnaavverr/STUDY/Research in Optimization/' ...
                '1. 0th order optimization/Figures_Sept/quartic'];
elseif  strcmp(ProbType,'RSB')
    filepath = ['E:/Google Drive_nnaavverr/STUDY/Research in Optimization/' ...
                '1. 0th order optimization/Figures_Sept/Rosenbrock'];
elseif  strcmp(ProbType,'SKQ')
    filepath = ['E:/Google Drive_nnaavverr/STUDY/Research in Optimization/' ...
                '1. 0th order optimization/Figures_Sept/skew_quartic'];
elseif  strcmp(ProbType,'HEQ')
    filepath = ['E:/Google Drive_nnaavverr/STUDY/Research in Optimization/' ...
                '1. 0th order optimization/Figures_Sept/Heat_EQ'];
end

formattype = 'epsc';
linewidth = 2; % 4; for small # exp
% #Queries/n  vs ftn val
fig = figure();
hold on;
if ~displog
    for i=1:N
        slen = length(Results{i}.name);
        if strcmp(Results{i}.name(end-1:end),'QN')
            plot(Results{i}.num_queries / n, (Results{i}.objval_seq), ':*', 'LineWidth',linewidth);
        elseif (slen >= 8 && strcmp(Results{i}.name(1:8),'Rand Inv'))
            plot(Results{i}.num_queries / n, (Results{i}.objval_seq), 'm', 'LineWidth',linewidth);
        elseif slen>=21 && strcmp(Results{i}.name(1:21),'Newton-type Orth Samp')
            plot(Results{i}.num_queries / n, (Results{i}.objval_seq), ':', 'LineWidth',linewidth);
        else
            plot(Results{i}.num_queries / n, (Results{i}.objval_seq), 'LineWidth',linewidth);
        end
    end
    ylabel('f(x)');
    title([fparam.name, ': (Number of Queries / dim) vs f(x), dim = ',num2str(n)]);
else
    for i=1:N
        slen = length(Results{i}.name);
        if isfield(Results{i}, 'rmls')
            if (Results{i}.rmls(2) == Results{i}.rmls(1)) && (Results{i}.rmls(3) == 1)
                plot(Results{i}.num_queries / n, log10(Results{i}.objval_seq - fmin), ':*', 'LineWidth',linewidth);
            elseif (Results{i}.rmls(2) == Results{i}.rmls(1)) && (Results{i}.rmls(3) > 1)
                plot(Results{i}.num_queries / n, log10(Results{i}.objval_seq - fmin), ':', 'LineWidth',linewidth);
            elseif (Results{i}.rmls(2) > Results{i}.rmls(1))
                plot(Results{i}.num_queries / n, log10(Results{i}.objval_seq - fmin), '-.', 'LineWidth',linewidth);
            else
                plot(Results{i}.num_queries / n, log10(Results{i}.objval_seq - fmin), '--', 'LineWidth',linewidth);
            end
        elseif strcmp(Results{i}.name(end-1:end),'QN')
            plot(Results{i}.num_queries / n, log10(Results{i}.objval_seq - fmin), ':*', 'LineWidth',linewidth);
        elseif (slen >= 8 && strcmp(Results{i}.name(1:8),'Rand Inv'))
            plot(Results{i}.num_queries / n, log10(Results{i}.objval_seq - fmin), 'm', 'LineWidth',linewidth);
        elseif slen>=21 && strcmp(Results{i}.name(1:21),'Newton-type Orth Samp')
            plot(Results{i}.num_queries / n, log10(Results{i}.objval_seq - fmin), ':', 'LineWidth',linewidth);
        else
            plot(Results{i}.num_queries / n, log10(Results{i}.objval_seq - fmin), 'LineWidth',linewidth);
        end
    end
     ylabel('log10(|f(x_k)-f_{min}|)');
     title([fparam.name, ': (Number of Queries / dim) vs log10(|f(x_k)-f_{min}|), dim = ',num2str(n)]);
end

legend(legend_str, 'location','best');
% ylabel('log10(|f(x_k)-f_{min}|)');
xlabel('Number of Queries / dim');
% savefig(fullfile(filepath,[ProbType,'_numq_vs_fval.fig']));
% saveas(fig, fullfile(filepath,[ProbType,'_numq_vs_fval']),formattype);
% saveas(fig, fullfile(filepath,[ProbType,'_numq_vs_fval']),'png');
% legend(legend_str, 'location','best','Orientation','horizontal');

% % #iterations vs ftn val
% fig = figure();
% hold on;
% if ~displog
%     for i=1:N
%         slen = length(Results{i}.name);
%         if strcmp(Results{i}.name(end-1:end),'QN')
%             plot((Results{i}.objval_seq), '-*', 'LineWidth',linewidth);
%         elseif (slen >= 8 && strcmp(Results{i}.name(1:8),'Rand Inv'))
%             plot((Results{i}.objval_seq), 'm', 'LineWidth',linewidth);
%         elseif slen>=21 && strcmp(Results{i}.name(1:21),'Newton-type Orth Samp')
%             plot((Results{i}.objval_seq ), ':', 'LineWidth',linewidth);
%         else
%             plot((Results{i}.objval_seq), 'LineWidth',linewidth);
%         end
%     end
%     ylabel('f(x)');
% else
%     for i=1:N
%         slen = length(Results{i}.name);
%         if strcmp(Results{i}.name(end-1:end),'QN')
%             plot(log10(Results{i}.objval_seq - fmin), '-*', 'LineWidth',linewidth);
%         elseif (slen >= 8 && strcmp(Results{i}.name(1:8),'Rand Inv'))
%             plot(log10(Results{i}.objval_seq - fmin), 'm', 'LineWidth',linewidth);
%         elseif slen>=21 && strcmp(Results{i}.name(1:21),'Newton-type Orth Samp')
%             plot(log10(Results{i}.objval_seq - fmin), ':', 'LineWidth',linewidth);
%         else
%             plot(log10(Results{i}.objval_seq - fmin), 'LineWidth',linewidth);
%         end
%     end
%     ylabel('log10(|f(x_k)-f_{min}|)');
% end
% title([fparam.name, ': Number of iterations vs log10(|f(x_k)-f_{min}|), dim = ',num2str(n)]);
% legend(legend_str, 'location','best');
% xlabel('Number of iterations');
% % savefig(fullfile(filepath,[ProbType,'_numit_vs_fval.fig']));
% % saveas(fig, fullfile(filepath,[ProbType,'_numit_vs_fval']),formattype);
% % saveas(fig, fullfile(filepath,[ProbType,'_numit_vs_fval']),'png');
% legend(legend_str, 'location','best');
% % legend(legend_str, 'location','best','Orientation','horizontal');

if isfield(fparam, 'fend')
    fig = figure();
    hold on;
    if ~displog
        for i=1:N
            slen = length(Results{i}.name);
            fmin = fparam.fend;
            fstart = Results{i}.objval_seq(1);
            if strcmp(Results{i}.name(end-1:end),'QN')
                plot((Results{i}.objval_seq - fmin)/(fstart-fmin), '-*', 'LineWidth',linewidth);
            elseif (slen >= 8 && strcmp(Results{i}.name(1:8),'Rand Inv'))
                plot((Results{i}.objval_seq - fmin)/(fstart-fmin), 'm', 'LineWidth',linewidth);
            elseif slen>=21 && strcmp(Results{i}.name(1:21),'Newton-type Orth Samp')
                plot((Results{i}.objval_seq - fmin)/(fstart-fmin), ':', 'LineWidth',linewidth);
            else
                plot((Results{i}.objval_seq - fmin)/(fstart-fmin), 'LineWidth',linewidth);
            end
        end
        ylabel('f(x)');
    else
        for i=1:N
            slen = length(Results{i}.name);
            fmin = fparam.fend;
            fstart = Results{i}.objval_seq(1);
            if strcmp(Results{i}.name(end-1:end),'QN')
                plot(Results{i}.num_queries, log10((Results{i}.objval_seq - fmin)/(fstart-fmin)), '-*', 'LineWidth',linewidth);
            elseif (slen >= 8 && strcmp(Results{i}.name(1:8),'Rand Inv'))
                plot(Results{i}.num_queries, log10((Results{i}.objval_seq - fmin)/(fstart-fmin)), 'm', 'LineWidth',linewidth);
            elseif slen>=21 && strcmp(Results{i}.name(1:21),'Newton-type Orth Samp')
                plot(Results{i}.num_queries, log10((Results{i}.objval_seq - fmin)/(fstart-fmin)), ':', 'LineWidth',linewidth);
            else
                plot(Results{i}.num_queries, log10((Results{i}.objval_seq - fmin)/(fstart-fmin)), 'LineWidth',linewidth);
            end
        end
        ylabel('log10(|f(x_k)-f_{min}|/|f(x_0)-f_{min}|)');
    end
    title([fparam.name, ': Number of queries vs log10(|f(x_k)-f_{min}|/|f(x_0)-f_{min}|), dim = ',num2str(n)]);
    legend(legend_str, 'location','best');
    xlabel('Number of iterations');
    % savefig(fullfile(filepath,[ProbType,'_numit_vs_fval.fig']));
    % saveas(fig, fullfile(filepath,[ProbType,'_numit_vs_fval']),formattype);
    % saveas(fig, fullfile(filepath,[ProbType,'_numit_vs_fval']),'png');
    legend(legend_str, 'location','best');
    % legend(legend_str, 'location','best','Orientation','horizontal');
end

% % #iterations vs |H - H_k| for QN only
% % 
% % if strcmp(Results{1}.name(end-1:end),'QN') && isfield(fparam,'Hess')
% %     fig = figure();
% %     hold on;
% %     Hrelnorms = zeros(length(Results{1}.Hseq),1);
% %     Hrelnorms2 = zeros(length(Results{end}.Hseq),1);
% %     for i=1:length(Results{1}.Hseq)
% %         Hrelnorms(i) = norm(Results{1}.Hseq{i} - fparam.Hess,'fro')/norm(fparam.Hess,'fro');
% %     end
% %     plot(Results{1}.num_queries(2:end) / n, Hrelnorms, 'LineWidth',linewidth);
% %     for i=1:length(Results{end}.Hseq)
% %         Hrelnorms2(i) = norm(Results{end}.Hseq{i} - fparam.Hess,'fro')/norm(fparam.Hess,'fro');
% %     end
% %     plot(Results{end}.num_queries(2:end) / n, Hrelnorms2, 'LineWidth',linewidth);
% %     title(['Number of queries vs  |H_k - H_{opt}|_F / |H_{opt}|_F, dim = ', num2str(n)]);
% %     legend('QN', 'Subspace QN', 'location','best','Orientation','horizontal');
% %     savefig(fullfile(filepath,[ProbType,'_numit_vs_HessErr.fig']));
% %     saveas(fig, fullfile(filepath,[ProbType,'_numit_vs_HessErr']),formattype);
% % end
% 
% 
% fig = figure();
% hold on;
% solerrs = cell(N,1);
% for j=1:N
%     slen = length(Results{j}.name);
%     l = length(Results{j}.sol);
%     solerrs{j} = zeros(l,1);
%     for i=1:l
%         solerrs{j}(i) = norm(Results{j}.sol{i}-true_sol); % norm(true_sol) = 1;
%     end
%     if strcmp(Results{j}.name(end-1:end),'QN') || (slen >= 8 && strcmp(Results{j}.name(1:8),'Rand Inv'))
%         plot(Results{j}.num_queries / n, log10(solerrs{j}), '-*', 'LineWidth',linewidth);
%     elseif slen>=21 && strcmp(Results{j}.name(1:21),'Newton-type Orth Samp')
%         plot(Results{j}.num_queries / n, log10(solerrs{j}), ':', 'LineWidth',linewidth);
%     else
%         plot(Results{j}.num_queries / n, log10(solerrs{j}), 'LineWidth',linewidth);
%     end
% end
% title(['Number of queries / dim vs  log(|x_k - x^*|), dim = ', num2str(n)]);
% legend(legend_str, 'location','best');
% savefig(fullfile(filepath,[ProbType,'_numq_vs_SolErr.fig']));
% saveas(fig, fullfile(filepath,[ProbType,'_numq_vs_SolErr']),formattype);
% saveas(fig, fullfile(filepath,[ProbType,'_numq_vs_SolErr']),'png');
% 
% % output
% dispres = struct;
% dispres.err = err;
% 
% %%
% %TEST
% 
% 
% fig = figure();
% hold on;
% for i=1:N
%     slen = length(Results{i}.name);
%     if strcmp(Results{i}.name(end-1:end),'QN') || (slen >= 8 && strcmp(Results{i}.name(1:8),'Rand Inv'))
%         plot(log10( 1 + Results{i}.num_queries / n), log10(Results{i}.objval_seq - fmin), '-*', 'LineWidth',linewidth);
%     elseif slen>=21 && strcmp(Results{i}.name(1:21),'Newton-type Orth Samp')
%          plot(log10( 1 + Results{i}.num_queries / n), log10(Results{i}.objval_seq - fmin), ':', 'LineWidth',linewidth);
%     else
%         plot(log10( 1 + Results{i}.num_queries / n), log10(Results{i}.objval_seq - fmin), 'LineWidth',linewidth);
%     end
% end
% title(['log10(1 + Number of Queries / dim) vs log10(|f(x_k)-f_{min}|), dim = ',num2str(n)]);

end