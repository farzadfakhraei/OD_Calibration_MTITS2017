function  [DATA_OUTPUT,DATA_INPUT,DATA_STRUCT]= part_B(dt,totT,ODmatrices,timeSeries,...
    links,nodes,det,DATA_INPUT)
if DATA_INPUT.partB.method==1
    approach = 'All OD pairs at once';
else
    approach = 'splitting into two clusters';
end

%%%%% note that %%%%%
% in this experiment we fist checked the approximated hessian to decide
% having two separate clusters although in the code we assumed this
% information is already know. as you run the code you will get the figure
% which shows approximated hessian
fprintf(1,'\n\nPart III.B of the paper with %s approach and Scenario %d\n',...
    approach, DATA_INPUT.partB.scenario)


%build the full ODmatrix
[ODmatrix,origins,destinations] = buildODmatrix(ODmatrices,timeSeries,dt,totT);

%% Computation of the base scenario
% First the base scenario is computed with a cold start. The iteration
% procedure of each time interval is initialized with CVN values of the
% previous time slice. This assumes that initially no vehicles proceed to
% the next time interval.
%
disp('Running I-LTM with a cold start')
[ODmatrix,DATA_STRUCT,DATA_OUTPUT] = Base_Senario(links,nodes,origins,destinations,ODmatrix,dt,totT,det);


%% Reducing Network dimension
% Here we only select two destinations, all origins going to destination 6
% and 16 and later analysis consider only these OD pairs
%
[odpairs_ind_des_based,~] = OD_pair_selection(sum(ODmatrix,3));
OD_selection{1} = [odpairs_ind_des_based{3}; odpairs_ind_des_based{8}];
DATA_OUTPUT.experiments = OD_selection;
DATA_OUTPUT.selectedExperiment = OD_selection{1};
O = OD_selection{1}(:,1)';
D =unique(OD_selection{1}(:,2))';
DATA_INPUT.origin_nodes_index = O;
DATA_INPUT.destination_nodes_index = D;
DATA_INPUT.origin_nodes = origins(1,O(1,1:end));
DATA_INPUT.destination_nodes = destinations(1,D(1,1:end));

%% Compute the response of a initial matrix perturbation 
disp('Running I-LTM with a warm start after purturbing the True matrix')
[odmatrix_perturbed,DATA_STRUCT,DATA_OUTPUT] = Perturbed_Senario(ODmatrix,DATA_INPUT,DATA_STRUCT,DATA_OUTPUT);

%% Visual definition of OD clusters 
% Here we use the perturbed matrix to define the clusters based on
% visual inspection of approximated Hessian 
%
[~,~,x_0,fixed,DATA_STRUCT]  = OD_initial_point_selection(odmatrix_perturbed,DATA_INPUT,DATA_STRUCT);
Xn{1}=x_0; Xn{2} = fixed;    
[exp_jac,~,~,deviation,DATA_STRUCT]= Calculate_Jacobian(Xn,DATA_STRUCT);
iterator = 1;
start_col = 1;
numberOfOD= size(exp_jac,2)/DATA_STRUCT.totT;
for end_col =numberOfOD:numberOfOD:size(exp_jac,2)
    temp_A = exp_jac(:,start_col:end_col);
        appro_hessian{iterator} = temp_A'*temp_A;
    start_col = start_col +numberOfOD;
    iterator = iterator +1;
end

%This was used initially as visual indication for the clustering
heatmap(appro_hessian{1})

%% initialization
DATA_STRUCT.iterations = 1;
DATA_OUTPUT.initials.seed_matrix = odmatrix_perturbed;
[~,DATA_OUTPUT.simFlow_perturb,DATA_STRUCT]= Objective_function(Xn,DATA_STRUCT);
DATA_OUTPUT.measurement_error(:,DATA_STRUCT.iterations) = rmse(reshape(DATA_STRUCT.Y_real,[],1),reshape(DATA_OUTPUT.simFlow_perturb,[],1));
RMSE_error = DATA_OUTPUT.measurement_error;
congestion_up{1} = DATA_STRUCT.con_up;
congestion_down{1} = DATA_STRUCT.con_down;

%Levenberg-Marquardt parameteres
vLM = 0.1;
OD_selection{1} = [odpairs_ind_des_based{3}];
OD_selection{2} = [odpairs_ind_des_based{8}];
global_termination = false;
internal_iteration = 0;
lambdaLM_int = ones(size(OD_selection,2),1)';

%Print the Initial Conditions
fprintf(1,'\n\nInitial Conditions are:\n')
fprintf(1,'Initial Lambda:\t%6f\n',lambdaLM_int)
fprintf(1,'Initial RMSE error:\t%0.2f\n',RMSE_error)
fprintf(1,'Problem dimention:\t%g\n',length(Xn{1}))
tStart = tic;

%define whether all OD pairs would be calibrated at once or two clusters
switch DATA_INPUT.partB.method
    case 1
        clusering = false;
        possible_clusters = 1;
    case 2
        clusering = true; 
        possible_clusters = 2;
end

%main algorithm
while (~global_termination)
    internal_iteration = internal_iteration+1;
    
    for zz = 1:possible_clusters
        
        inner_iteration = 0;
        
        if clusering
            % initialize lambda
            lambdaLM = lambdaLM_int(zz) ;
            
            % in the case of clustering, after each convergence we switch
            % between each cluster
            OD_selection{zz}  = OD_selection{zz};
            O = [OD_selection{zz}(:,1)]';
            D =unique([OD_selection{zz}(:,2)]);
            DATA_INPUT.origin_nodes_index = O;
            DATA_INPUT.destination_nodes_index = D;
            DATA_INPUT.origin_nodes = origins(1,O(1,1:end));
            DATA_INPUT.destination_nodes = destinations(1,D);
            
            if (mod(zz,2) ==0 || internal_iteration>1)
                odmatrix_perturbed = reshape(DATA_OUTPUT.inter_ODmatrix{1,end},size(origins,2),size(destinations,2),totT);
            end
            
            %Define the initial points
            [~,~,x_0,fixed,DATA_STRUCT]  = OD_initial_point_selection(odmatrix_perturbed,DATA_INPUT,DATA_STRUCT);
            Xn{1}=x_0; Xn{2} = fixed;       
            DATA_INPUT.origin_nodes_index = O;
            DATA_INPUT.destination_nodes_index = D;
            DATA_INPUT.origin_nodes = origins(1,O(1,1:end));
            DATA_INPUT.destination_nodes = destinations(1,D(1,1:end));
        else
            lambdaLM = lambdaLM_int(zz) ;
        end
        converged = false;
        
        % main Loop
        while (~converged)
            DATA_STRUCT.iterations = DATA_STRUCT.iterations +1;
            inner_iteration = inner_iteration + 1;
            
            %calculate of the jacobian
            if DATA_STRUCT.iterations==2
                [Jac,~,~,deviation,DATA_STRUCT]= Calculate_Jacobian(Xn,DATA_STRUCT);
            else
                congestion_check = sum(sum(abs(congestion_up{DATA_STRUCT.iterations-1}- congestion_up{DATA_STRUCT.iterations-2})));
                if any(congestion_check)
                    
                    %calculating jacobian
                    [Jac,~,~,deviation,DATA_STRUCT]= Calculate_Jacobian(Xn,DATA_STRUCT);
                else
                    
                    %restore the jacobian
                    Jac = DATA_OUTPUT.jacobian{DATA_STRUCT.iterations-1} ;
                end
            end
            DATA_OUTPUT.jacobian{DATA_STRUCT.iterations} = sparse(Jac);
            
            
            %calculation of the search direction
            step_m = mldivide(transpose(Jac)*Jac + lambdaLM * diag(diag(transpose(Jac)*Jac)),(transpose(Jac)*deviation.measurement));  %normal LM
            DATA_OUTPUT.step{DATA_STRUCT.iterations} = step_m;
            
            %stop if wrong step is calculated 
            if any(isnan(step_m))
                step_m = zeros(length(step_m),1);
                converged = true;
            end
            
            %update to the new point
            temp_Xn = Xn{1};
            Xn{1} = Xn{1}+ step_m; 
            
            %hardcore negative values removal
            Xn{1}(Xn{1}<0) =1;
            
            %evaluate the new step
            [r,flow,DATA_STRUCT]= Objective_function(Xn,DATA_STRUCT);
            
            %evaluate the OF
            RMSE_error = rmse(reshape(DATA_STRUCT.Y_real,[],1),flow);

            %improvment check
            if RMSE_error< DATA_OUTPUT.measurement_error(:,DATA_STRUCT.iterations-1)
                
                f = 0.5 * (r.measurement' * r.measurement);
                g = Jac' * r.measurement;
                f_norm = norm(f);
                g_norm = norm(g);
                step_norm = norm(step_m);
                
                %display the results
                fprintf(1,'\n\nIteration:%g\n',DATA_STRUCT.iterations) %#ok<*PRTCAL>
                fprintf(1,'Step Norm:\t%0.5f\n',step_norm)
                fprintf(1,'Residual Norm:\t%0.5f\n',f_norm)
                fprintf(1,'Norm of Gradient:\t%0.5f\n',g_norm)
                fprintf(1,'Lambda:\t%4e\n',lambdaLM)
                fprintf(1,'RMSE error is:\t%0.5f\n',RMSE_error)
                fprintf(1,'Problem dimention:\t%g\n',length(Xn{1}))
                
%                 check for the damp factor update KWAK
%                 implementation
                if DATA_STRUCT.iterations>2
                    previous_step = DATA_OUTPUT.step{:,DATA_STRUCT.iterations-1};
                    wip = (dot(step_m,previous_step)) /(norm(step_m)*norm(previous_step));
                    lambdaLM = lambdaLM * (vLM^wip);
                end
            else
                % IF OF increases, we need to take smaller step. 
                internalcount = 1;
                fprintf('\n\n Wrong direction, Lambda would increase\n')
                while RMSE_error > DATA_OUTPUT.measurement_error(:,DATA_STRUCT.iterations-1) && internalcount< 10
                    %If the error has increased as a result the update, then retract the step (i.e. reset the weights to
                    % their previous values) and increase ? by a factor of 10 or some such significant factor.
                    Definite = transpose(Jac)*Jac + lambdaLM *  diag(diag(transpose(Jac)*Jac));
                    [~,positive_definite] = chol(Definite);
                    if positive_definite == 0
                        lambdaLM = lambdaLM /vLM;
                    else
                        for r = 1:size(Definite,1)
                            d(r,:)= sum(abs(Definite(r,:)));
                        end
                        lambdaLM = min(d);
                    end
                    if lambdaLM>1.0000e+3
                        lambdaLM = 1.0000e+3;
                    end
                    %calculate the direction
                    step_m = mldivide(transpose(Jac)*Jac + lambdaLM *  diag(diag(transpose(Jac)*Jac)) ,(transpose(Jac)*deviation.measurement));  %normal LM
                    
                    %update to the new point
                    Xn{1} = temp_Xn;
                    Xn{1} = Xn{1}+ step_m;
                    
                    %hardcore negative values removal
                    Xn{1}(Xn{1}<0) =1;
                    
                    %evaluate the new step
                    [~,flow,DATA_STRUCT]= Objective_function(Xn,DATA_STRUCT);
                    
                    %display the error on detectors
                    RMSE_error = rmse(reshape(DATA_STRUCT.Y_real,[],1),flow);
                    fprintf(1,'Number of iteration in the wrong step:%g\n',internalcount)
                    internalcount = internalcount +1;
                    
                end
                fprintf(1,'Lambda reduced, RMSE error is:\t%0.2f\n',RMSE_error)
                fprintf(1,'New Lambda value is:\t%4e\n',lambdaLM)
                DATA_STRUCT.iterations = DATA_STRUCT.iterations - 1;
            end
            
            %% save intermediate results
            DATA_OUTPUT.iterations = DATA_STRUCT.iterations;
            [DATA_OUTPUT]= store_intermediate_results(Xn,RMSE_error,DATA_STRUCT,DATA_OUTPUT,DATA_INPUT);
            error_matrix = rmse(reshape(DATA_STRUCT.odmatrix_original,[],1),reshape(DATA_OUTPUT.inter_ODmatrix{DATA_STRUCT.iterations},[],1));   %showing current RMSE error
            display(['OD Matrix error: ', num2str(error_matrix)]);
            congestion_up{DATA_STRUCT.iterations} = DATA_STRUCT.con_up;
            congestion_down{DATA_STRUCT.iterations} = DATA_STRUCT.con_down;
            
            %% Check for convergence
            stopping_criteria = (((DATA_OUTPUT.measurement_error(:,DATA_STRUCT.iterations-1) - ...
                RMSE_error)/DATA_OUTPUT.measurement_error(:,DATA_STRUCT.iterations-1)))*100;
            fprintf(1,'Percentage improvment rate:\t%0.5f\n',stopping_criteria)
            tElapsed = toc(tStart);
            if  (stopping_criteria < 1 || tElapsed > 200)
                converged = true;
            end
            
            
        end
        lambdaLM_int(zz)= lambdaLM ;
        fprintf('One cluster is done, go to the next one\n');
        tElapsed = toc(tStart);
        DATA_OUTPUT.total_time_spend_sec = tElapsed;
        fprintf(1,'\n Total Time Spent(sec): %g\n',tElapsed)
        
    end
    
    %stop if budget time exceed 200 seconds
    if (tElapsed > 200)
        global_termination = true;
    end
end
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%            Functions                %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation of the base scenario
function [ODmatrix,DATA_STRUCT,DATA_OUTPUT] = Base_Senario(links,nodes,origins,destinations,ODmatrix,dt,totT,det)
% First the base scenario is computed with a cold start. The iteration
% procedure of each time interval is initialized with CVN values of the
% previous time slice. This assumes that initially no vehicles proceed to
% the next time interval.

%set faster lookup structures
[links,node_prop] = dataParser(links,nodes,origins,destinations,dt);

%set turning fractions faster (based on free flow conditions)
TF=TF_init(node_prop,links,destinations,dt,totT);

%run I-LTM
tic
[cvn_up_t,cvn_down_t,con_up_t,con_down_t] = ILTM_cold_mex(node_prop,links,origins,destinations,ODmatrix,dt,totT,TF);
toc

%compute the density
cvn_up_true=reshape(sum(cvn_up_t,2),[],totT+1);
cvn_down_true=reshape(sum(cvn_down_t,2),[],totT+1);
[simDensity_true] = cvn2dens(cvn_up_true,cvn_down_true,totT,links);
simFlow_true_all = cvn2flows(cvn_down_true,dt);
[simFlow_true] = simFlow_true_all(det,:);
[cong_true] = con_down_t(det,:);

%structure of OD pairs definition
temp_OD_mat = sum(ODmatrix,3);
nod=nnz(temp_OD_mat);
[I,J] = ind2sub(size(temp_OD_mat),find(temp_OD_mat));
odpairs_ind=[I,J];
sp_i=repmat(odpairs_ind(:,1),totT,1);
sp_j=repmat(odpairs_ind(:,2),totT,1);
sp_k=sort(reshape([1:totT]'*ones(1,nod),1,[])');

odmatrix_original = ODmatrix;
DATA_STRUCT.prev_x_with_zeros = reshape(odmatrix_original,1,[])';
DATA_STRUCT.cvn_up = cvn_up_t;
DATA_STRUCT.cvn_down = cvn_down_t;
DATA_STRUCT.con_up = con_up_t;
DATA_STRUCT.con_down = con_down_t;

DATA_STRUCT.sp_i = sp_i;
DATA_STRUCT.sp_j = sp_j;
DATA_STRUCT.sp_k = sp_k;

DATA_STRUCT.node_prop = node_prop;
DATA_STRUCT.nodes = nodes;
DATA_STRUCT.links = links;
DATA_STRUCT.origins = origins;
DATA_STRUCT.destinations = destinations;

DATA_STRUCT.dt = dt;
DATA_STRUCT.totT = totT;
DATA_STRUCT.TF = TF;

DATA_STRUCT.det = det;
DATA_STRUCT.Y_real = simFlow_true;
DATA_STRUCT.Y_con = con_down_t;
DATA_STRUCT.history = [];
DATA_OUTPUT.simFlow_true_all = simFlow_true_all;
DATA_OUTPUT.simFlow_true = simFlow_true;
DATA_OUTPUT.cong_true=cong_true;
DATA_OUTPUT.simDensity_true =simDensity_true;
DATA_STRUCT.odmatrix_original= odmatrix_original;

end

%% OD pairs Selection based on shared Destination or two by two
function [odpairs_ind_des_based,experiment] = OD_pair_selection(temp_OD_mat)
%extract OD pair information and put them in structure form
experiment5_OD_set = [];
m =1 ;
while m<=size(temp_OD_mat,2)
    temp_OD_mat_destination_based = temp_OD_mat(:,m);
    [I1,J1] = ind2sub(size(temp_OD_mat_destination_based),find(temp_OD_mat_destination_based));
    odpairs_ind_des_based{m}=[I1,m*J1];
    m = m+1;
end

%define the experiment tests where two origins share one destination
m = 1;
i = 1;
while m<= size(odpairs_ind_des_based,2)
    if size(odpairs_ind_des_based{1,m},1)==2
        temp_experi{i} = odpairs_ind_des_based{1,m};
        i = i +1;
    else
        for j =1:size(odpairs_ind_des_based{1,m},1)-1
            temp_experi{i} = odpairs_ind_des_based{1,m}(j:j+1,:);
            i = i +1;
        end
    end
    m = m+1;
end

%this section will be used after one time running all the experiments
for i = 1:length(temp_experi)
    temp_experi_reverse{1,i}= temp_experi{1,i}([2 1],:);
end
experiment = [temp_experi,temp_experi_reverse];
end

%% define the starting point of OD matrix variables
function [ub,lb,x_0,x_0_fix,DATA_STRUCT] = OD_initial_point_selection(odmatrix_s,DATA_INPUT,DATA_STRUCT)
O = DATA_INPUT.origin_nodes_index;
D = DATA_INPUT.destination_nodes_index;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed OD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the ones we want to keep fixed
temp_OD_fixed = odmatrix_s;
temp_OD_fixed(temp_OD_fixed>0)=1;
for i = 1:length(DATA_INPUT.origin_nodes)
    for j = 1:length(DATA_INPUT.destination_nodes)
        if (odmatrix_s(O(i),D(j),:))>0
            temp_OD_fixed(O(i),D(j),:) =0;
        end
    end
end

%define the indexes of the OD pairs we want to optimize
freezed_OD = reshape(temp_OD_fixed,length(DATA_STRUCT.origins),length(DATA_STRUCT.destinations),DATA_STRUCT.totT');
temp_OD_mat = sum(freezed_OD,3);
nod=nnz(temp_OD_mat);
[I,J] = ind2sub(size(temp_OD_mat),find(temp_OD_mat));
odpairs_ind=[I,J];
DATA_STRUCT.sp_i_f=repmat(odpairs_ind(:,1),DATA_STRUCT.totT',1);
DATA_STRUCT.sp_j_f=repmat(odpairs_ind(:,2),DATA_STRUCT.totT',1);
DATA_STRUCT.sp_k_f=sort(reshape([1:DATA_STRUCT.totT']'*ones(1,nod),1,[])');

%restore the OD matrix
OD_freezed = odmatrix_s.*temp_OD_fixed;
x_0_fix = reshape(OD_freezed,1,[])';

% %make initial values without zeros
x_0_fix(x_0_fix==0)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variable OD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the ones we want to optimize
temp_OD_optimize = 0*odmatrix_s;
for i = 1:length(DATA_INPUT.origin_nodes)
    for j = 1:length(DATA_INPUT.destination_nodes)
        if (odmatrix_s(O(i),D(j),:))>0
            temp_OD_optimize(O(i),D(j),:) =1;
        end
    end
end
%define the fixed OD
variable_OD = reshape(temp_OD_optimize,length(DATA_STRUCT.origins),length(DATA_STRUCT.destinations),DATA_STRUCT.totT');
temp_OD_mat = sum(variable_OD,3);
nod=nnz(temp_OD_mat);
[I,J] = ind2sub(size(temp_OD_mat),find(temp_OD_mat));
odpairs_ind=[I,J];
% odpairs_ind=[O,D];
DATA_STRUCT.sp_i_v=repmat(odpairs_ind(:,1),DATA_STRUCT.totT',1);
DATA_STRUCT.sp_j_v=repmat(odpairs_ind(:,2),DATA_STRUCT.totT',1);
DATA_STRUCT.sp_k_v=sort(reshape([1:DATA_STRUCT.totT']'*ones(1,nod),1,[])');

%restore the OD matrix
OD_variable = odmatrix_s.*temp_OD_optimize;
x_0 = reshape(OD_variable,1,[])';
x_0(x_0==0)=[];
lb = zeros(size(x_0))+.001;
ub = x_0 *5;
ub(ub==0)=[];

end

%% Compute the response of a initial matrix purtubation (if Ture matrix is available)
function [odmatrix_perturbed,DATA_STRUCT,DATA_OUTPUT] = Perturbed_Senario(ODmatrix,DATA_INPUT,DATA_STRUCT,DATA_OUTPUT)
% Here the base original OD matrix is purturb based on different
% purtubation criteria, then the purturbed matrix and simulation results
% are stored for later comparision


odmatrix_perturbed = ODmatrix;

switch DATA_INPUT.partB.scenario
    case 1
        odmatrix_perturbed(DATA_INPUT.origin_nodes_index(:),DATA_INPUT.destination_nodes_index(1),:)...
            = ((0.3)*rand(1))* ODmatrix(DATA_INPUT.origin_nodes_index(:),DATA_INPUT.destination_nodes_index(1),:);
        
        odmatrix_perturbed(DATA_INPUT.origin_nodes_index(:),DATA_INPUT.destination_nodes_index(2),:)...
            = ((0.7)*rand(1))* ODmatrix(DATA_INPUT.origin_nodes_index(:),DATA_INPUT.destination_nodes_index(2),:);
    case 2
        %OD pairs passing bottelneck starts far from solution
        odmatrix_perturbed(DATA_INPUT.origin_nodes_index(:),DATA_INPUT.destination_nodes_index(1),:)...
            = ((0.7)*rand(1))* ODmatrix(DATA_INPUT.origin_nodes_index(:),DATA_INPUT.destination_nodes_index(1),:);
        
        odmatrix_perturbed(DATA_INPUT.origin_nodes_index(:),DATA_INPUT.destination_nodes_index(2),:)...
            = ((0.3)*rand(1))* ODmatrix(DATA_INPUT.origin_nodes_index(:),DATA_INPUT.destination_nodes_index(2),:);
end

x_with_zeros=reshape(odmatrix_perturbed,1,[])';

%check if any OD need to be updated
if any(DATA_STRUCT.prev_x_with_zeros~=x_with_zeros)
    [i,~,k]=ind2sub(size(odmatrix_perturbed),find(DATA_STRUCT.prev_x_with_zeros~=x_with_zeros));
    nodes2update = false(length(DATA_STRUCT.nodes.id),DATA_STRUCT.totT+1);
    for i_ind=i'
        for k_ind=k'
            nodes2update(DATA_STRUCT.origins(i_ind),k_ind+1) = true;
        end
    end
    [cvn_up,cvn_down,con_up,con_down]=ILTM_warm_mex(DATA_STRUCT.node_prop,DATA_STRUCT.links,DATA_STRUCT.origins,DATA_STRUCT.destinations,odmatrix_perturbed,DATA_STRUCT.dt,DATA_STRUCT.totT,DATA_STRUCT.TF,DATA_STRUCT.cvn_up,DATA_STRUCT.cvn_down,DATA_STRUCT.con_up,DATA_STRUCT.con_down,nodes2update);
    
    %store the output of the simluation
    DATA_STRUCT.cvn_up = cvn_up;
    DATA_STRUCT.cvn_down = cvn_down;
    DATA_STRUCT.con_up = con_up;
    DATA_STRUCT.con_down = con_down;
    DATA_STRUCT.prev_x_with_zeros = x_with_zeros;
end


cvn_up_perturb=reshape(sum(DATA_STRUCT.cvn_up,2),[],DATA_STRUCT.totT+1);
cvn_down_perturb=reshape(sum(DATA_STRUCT.cvn_down,2),[],DATA_STRUCT.totT+1);
simDensity_perturbed = cvn2dens(cvn_up_perturb,cvn_down_perturb,DATA_STRUCT.totT,DATA_STRUCT.links);
simFlow_perturbed_all = cvn2flows(cvn_down_perturb,DATA_STRUCT.dt);
simFlow_perturbed = simFlow_perturbed_all(DATA_STRUCT.det,:);
cong_perturbed = DATA_STRUCT.con_down(DATA_STRUCT.det,:);

DATA_STRUCT.odmatrix_purturbed = odmatrix_perturbed;
DATA_OUTPUT.simFlow_perturb_all = simFlow_perturbed_all;
DATA_OUTPUT.simFlow_perturb = simFlow_perturbed;
DATA_OUTPUT.cong_perturbed=cong_perturbed;
DATA_OUTPUT.simDensity_perturb =simDensity_perturbed;
end

%% calculate the jacobian
function [jac,grads,flow_current,deviation,DATA_STRUCT] = Calculate_Jacobian(Xn,DATA_STRUCT)

if size(Xn,2)>1
    x = Xn{1};
    fixed = Xn{2};
else
    x = Xn;
end
[deviation,flow_current,DATA_STRUCT]=Objective_function(Xn,DATA_STRUCT);                % function value
n=numel(x);                     % size of independent
m=numel(deviation.measurement);      % size of dependent

jac=zeros(m,n);                   % allocate memory for the Jacobian matrix(initial usaged)
h=1;                            % differentiation step size (vehcile)
parfor i=1:n                       % loop for each independent variable
    xf=x;                       % reference point
    xf(i)=xf(i)+h;              % calculate the new point forward
    if size(Xn,2)>1
        Xnf{1}= xf; Xnf{2}= fixed;
    else
        Xnf= xf;
    end
    [~,Yf,~] = Objective_function(Xnf,DATA_STRUCT);    % calculate the next point
    jac(:,i)=(Yf-flow_current)/h;            % step differentiation forward finite differences (initial usaged)
end

jac(find(abs(jac)<0.01))=0;
transf_matrix{1} = speye(size(x(:,1),1));


for i=1:size(x,1)
    grad_meas(i,1) = deviation.measurement'*(-jac(:,i));
    grad_target(i,1) = deviation.target_matrix'*(transf_matrix{1}(:,i));
    
end
grads.meas{1} = grad_meas(:,1);
grads.target{1} = grad_target(:,1);
grads.total{1}= [grad_meas+grad_target];
end

%% Objective function
function [deviation,f,DATA_STRUCT]= Objective_function(Xn,DATA_STRUCT)
if size(Xn,2)>1
    x = Xn{1};
    fixed = Xn{2};
else
    x = Xn;
end


%initialize the OD matrix
odmatrix=zeros(length(DATA_STRUCT.origins),length(DATA_STRUCT.destinations),DATA_STRUCT.totT);

%check the inputs for the function
if size(Xn,2) == 1 || (isempty(Xn{2}))
    odmatrix(sub2ind(size(odmatrix),DATA_STRUCT.sp_i,DATA_STRUCT.sp_j,DATA_STRUCT.sp_k))=x;
else
    odmatrix(sub2ind(size(odmatrix),DATA_STRUCT.sp_i_v,DATA_STRUCT.sp_j_v,DATA_STRUCT.sp_k_v))=x;
    odmatrix(sub2ind(size(odmatrix),DATA_STRUCT.sp_i_f,DATA_STRUCT.sp_j_f,DATA_STRUCT.sp_k_f))=fixed;
end
x_with_zeros=reshape(odmatrix,1,[])';

%check which OD node should be updated
if any(DATA_STRUCT.prev_x_with_zeros~=x_with_zeros)
    [i,~,k]=ind2sub(size(odmatrix),find(DATA_STRUCT.prev_x_with_zeros~=x_with_zeros));
    nodes2update = false(length(DATA_STRUCT.nodes.id),DATA_STRUCT.totT+1);
    for i_ind=i'
        for k_ind=k'
            nodes2update(DATA_STRUCT.origins(i_ind),k_ind+1) = true;
        end
    end
    [cvn_up,cvn_down,con_up,con_down]=ILTM_warm_mex(DATA_STRUCT.node_prop,DATA_STRUCT.links,DATA_STRUCT.origins,DATA_STRUCT.destinations,odmatrix,DATA_STRUCT.dt,DATA_STRUCT.totT,DATA_STRUCT.TF,DATA_STRUCT.cvn_up,DATA_STRUCT.cvn_down,DATA_STRUCT.con_up,DATA_STRUCT.con_down,nodes2update);
    
    %store the output of the simluation
    DATA_STRUCT.cvn_up = cvn_up;
    DATA_STRUCT.cvn_down = cvn_down;
    DATA_STRUCT.con_up = con_up;
    DATA_STRUCT.con_down = con_down;
    DATA_STRUCT.prev_x_with_zeros = x_with_zeros;
end

% Evaluate the goal function
cvn_down_p=reshape(sum(DATA_STRUCT.cvn_down,2),[],DATA_STRUCT.totT+1);
flow = cvn2flows(cvn_down_p,DATA_STRUCT.dt);
flow= flow(DATA_STRUCT.det,:);
f=reshape(flow,[],1);

current_odmatrix=zeros(length(DATA_STRUCT.origins),length(DATA_STRUCT.destinations),DATA_STRUCT.totT);
if size(Xn,2) == 1 || (isempty(Xn{2}))
    x_without_zeros = x_with_zeros; x_without_zeros(x_without_zeros==0)=[];
    current_odmatrix(sub2ind(size(current_odmatrix),DATA_STRUCT.sp_i,DATA_STRUCT.sp_j,DATA_STRUCT.sp_k))=x_without_zeros;
    current_od_vector = reshape(current_odmatrix,[],1);
    current_od_vector_without_zero = current_od_vector;  current_od_vector_without_zero(current_od_vector_without_zero==0)=[];
    original_od_vector= reshape(DATA_STRUCT.odmatrix_original,[],1); original_od_vector(original_od_vector==0)=[];
    original_od_vector_without_zero = original_od_vector; original_od_vector_without_zero(original_od_vector_without_zero==0)=[];
    
else
    current_odmatrix(sub2ind(size(current_odmatrix),DATA_STRUCT.sp_i_v,DATA_STRUCT.sp_j_v,DATA_STRUCT.sp_k_v))=x;
    %         current_odmatrix(sub2ind(size(current_odmatrix),DATA_STRUCT.sp_i_f,DATA_STRUCT.sp_j_f,DATA_STRUCT.sp_k_f))=fixed;
    current_od_vector = reshape(current_odmatrix,[],1);
    current_od_vector_without_zero = current_od_vector;  current_od_vector_without_zero(current_od_vector_without_zero==0)=[];
    variable_od_index = current_od_vector; variable_od_index(variable_od_index>0)=1;
    
    original_od_vector= reshape(DATA_STRUCT.odmatrix_original,[],1);
    original_od_variable_vector = original_od_vector.*variable_od_index;
    original_od_vector_without_zero = original_od_variable_vector; original_od_vector_without_zero(original_od_vector_without_zero==0)=[];
    
end

deviation.measurement = reshape(DATA_STRUCT.Y_real - flow,[],1);
deviation.target_matrix = original_od_vector_without_zero - current_od_vector_without_zero;
end

%% save intermediate results
function [DATA_OUTPUT] = store_intermediate_results(Xn,RMSE_error,DATA_STRUCT,DATA_OUTPUT,DATA_INPUT)
O = DATA_INPUT.origin_nodes_index;
D = DATA_INPUT.destination_nodes_index;
%save results
DATA_OUTPUT.intermediate_x(:,DATA_OUTPUT.iterations) = {Xn{1}};
DATA_OUTPUT.measurement_error(:,DATA_OUTPUT.iterations) = RMSE_error;

%save intermediate matrix
od_mat=zeros(length(DATA_STRUCT.origins),length(DATA_STRUCT.destinations),DATA_STRUCT.totT);

%check the inputs for the function
if size(Xn,2) == 1 || (isempty(Xn{2}))
    od_mat(sub2ind(size(od_mat),DATA_STRUCT.sp_i,DATA_STRUCT.sp_j,DATA_STRUCT.sp_k))=Xn{1};
else
    od_mat(sub2ind(size(od_mat),DATA_STRUCT.sp_i_f,DATA_STRUCT.sp_j_f,DATA_STRUCT.sp_k_f)) = Xn{2};
    od_mat(sub2ind(size(od_mat),DATA_STRUCT.sp_i_v,DATA_STRUCT.sp_j_v,DATA_STRUCT.sp_k_v)) = Xn{1};
end
DATA_OUTPUT.inter_ODmatrix{DATA_OUTPUT.iterations} = od_mat;
end

