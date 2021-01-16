function [Real_Payoffs,Sim_Payoffs,Mu,GBM_Vol,GARCH_Param, results] = American_Option_Pricing(numBasis,Tenors,N_Training_Paths,N_Test_Paths,Realized_Vol_Window,r,m)
% DEFINITIONS OF THE INPUT VARIABLES:
% 1) Tenors - tenors of the studied American options (in time units stored in Settings.Time_Unit).
% 2) Training_Fraction - fraction of the data to be used for reinforcement learning.
%                        Chronologically, the training segment precedes the test segment.
% 3) N_Training_Paths - number of training paths.
% 4) N_Test_Paths - number of test paths.
% 5) Realized_Vol_Window - the window corresponding to 1/2 of the exponential decay if using exponentially weighted realized
%                          volatility as a state component. When Realized_Vol_Window <= 1, realized volatility is not used.
% 6) r - riskless interest rate corresponding
% 7) m - specifies the algorithm (1: LSPI, 2: FQI, 3: LSM)

Settings.Path           = 'C:\Users\Frosya\Desktop\Mishal Thesis\Code';
Settings.Recent_to_Old  = 1;
Settings.Interpolate    = 1;
Settings.Log_Return_Cap = 0.08;
Settings.Time_Unit      = 'day';
Settings.dt             = 1 / 252; % Time unit in years. Must correspond to Settings.Time_Unit.
Settings.Min_Iter       = 3; 
Settings.Max_Iter       = 15; 
Settings.price_Tol      = 1e-2;
Settings.Det_Tol        = 1e-12;
Settings.Vol_Basis      = 1; % if 1, then use simple vol. If 0, use complex vol.
Settings.Preset_Ini_Vol = 0;
Settings.Methods        = {'LSPI', 'FQI', 'LSM'};
Settings.Data_Sources   = {'GBM', 'GARCH', 'Data'};
Settings.Warnings       = 'on';
Settings.numTrainTraj_realData = [600, 500, 500];
Settings.numTestTraj_realData  = [500, 450, 250];
Settings.numTrainDataPts_GBM   = [662, 625, 751];
Settings.Weight_Tol      = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOADING THE DATA & INITIALIZATION

StringData = textread('Registry.csv','%s','delimiter','\n'); % '%s' means string
Tickers = {};
for t = 1 : length(StringData)
    Tickers = [Tickers, StringData(t)]; % names of the stocks
end
Results_File        = [Settings.Path, sprintf('Results_%i', abs(Realized_Vol_Window))];
I                   = length(Tickers); % Number of underlying assets (29)
TEN                 = length(Tenors); % Number of tenors (3)
Max_Tenor           = max(Tenors);
Mu                  = zeros(I,1) * NaN;
GARCH_Param         = zeros(3,I) * NaN; % 1) Constant; 2) GARCH(1); 3) ARCH(1);
Spec                = garch(1,1);
% store results into variables
GBM_Vol             = zeros(TEN,I)   * NaN;
Real_Payoffs        = zeros(TEN,I,3) * NaN; % 3 is for Weights trained with 1) Real data 2) GBM 3) GARCH
Real_Payoffs_NonStd = zeros(TEN,I,3) * NaN; % transformed back into stock's original price scale
Sim_Payoffs         = zeros(TEN,I,2) * NaN; % 2 is for Weights trained with 1) GBM 2) GARCH
Sim_Payoffs_NonStd  = zeros(TEN,I,2) * NaN;
Real_Boundaries     = cell(TEN,I,3);
Real_Boundaries_Exh = cell(TEN,I,3);
Sim_Boundaries      = cell(TEN,I,2);
Sim_Boundaries_exh  = cell(TEN,I,2);
Weights             = cell(TEN,I,3);
Prices_Train        = zeros(TEN,I,3) * NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS OF EACH ASSET
for i = 13 % 1 : I % looping over all underlying stocks
    fprintf( '\n Underlying %s:', Tickers{i} );
    Underlying = xlsread([Tickers{i}, '.csv']);
    Underlying = Underlying(:,4); % only store adjusted closing prices (default order of Recent - Old)
    T = length(Underlying);       % # data points (i.e. # all day by day closing prices)
    if Settings.Recent_to_Old     % re-order data old to recent
        Underlying = Underlying(T:-1:1);
    end
    if Settings.Interpolate       % deals with missing data
        Underlying = Data_Interpolation(Underlying);
    end
    % figure,plot(Underlying,'x-','markerSize',10); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight;
    K = Underlying(1); % Non-standardised strike
    
    % ----- LOOP OVER 3 TENORS
    for ten = 2% 1 : TEN
        % ----- Fitting GBM and GARCH models to real training data.
        fprintf( '\n Estimating GBM from real data \n' );
        % GBM MODEL PARAM EXTRACTION
        T_Train = Settings.numTrainDataPts_GBM(ten); % Using first X num data pts for trainig models (diff for diff tenors)
        % --- log return way
        % calculating log of return log(price_t+1 / price_t)
        Real_Log_Return     = log(Underlying(2:T_Train)) - log(Underlying(1:(T_Train-1))); 
        Mu(i,1)             = mean(Real_Log_Return) / Settings.dt;      % drift
        GBM_Vol(ten,i)      = std(Real_Log_Return) / sqrt(Settings.dt); % volatility, sigma for GBM;
        Cent_Log_Return     = Real_Log_Return - Mu(i,1) * Settings.dt; % standardising
        if strcmpi(Tickers{i},'PFE')
            Cent_Log_Return = max(Cent_Log_Return, -Settings.Log_Return_Cap); % a cap on 'PFE' stock's return
        end
        % GARCH MODEL PARAM EXTRACTION
        [Coeff]             = estimate(Spec,Cent_Log_Return); 
        GARCH_Param(1,i)    = Coeff.Constant; % fills first row
        GARCH_Param(2,i)    = Coeff.GARCH{1};
        GARCH_Param(3,i)    = Coeff.ARCH{1};
        % --- GBM (fix seed since 3 algos meant to be trained/tested same trajs --> gives fair comparison of algos)
        GBM_Train           = Simulator(N_Training_Paths, Tenors(ten), 'GBM',   GBM_Vol(ten,i),   r, -1, Settings); % N_Training_Paths
        OutputName = 'TrajsGBM_Train';
        save(['results/' OutputName '.mat'], 'GBM_Train'); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight;
        % figure,plot(GBM_Train','x-'); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight;
        % title('GBM Data Train')
        GBM_Test            = Simulator(N_Test_Paths,     Tenors(ten), 'GBM',   GBM_Vol(ten,i),   r, -1, Settings); % N_Test_Paths
        % figure,plot(GBM_Test','x-'); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight;
        % title('GBM Data Test')
        % --- GARCH
        GARCH_Train         = Simulator(N_Training_Paths, Tenors(ten), 'GARCH', GARCH_Param(:,i), r, -1, Settings); 
        % figure,plot(GARCH_Train','x-'); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight;
        % title('GARCH Data Train')
        GARCH_Test          = Simulator(N_Test_Paths,     Tenors(ten), 'GARCH', GARCH_Param(:,i), r, -1, Settings); 
        % figure,plot(GARCH_Test','x-'); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight;
        % title('GARCH Data Test')
        % --- Set initial GARCH vols for the testing stage
        GARCH_Mean_Vol      = sqrt( GARCH_Param(1,i) / (1 - GARCH_Param(2,i) - GARCH_Param(3,i)) ) / sqrt(Settings.dt); % Note that the volatility is annualized.
        if Settings.Preset_Ini_Vol
            Ini_GBM_Vol     = GBM_Vol(ten,i);
            Ini_GARCH_Vol   = GARCH_Mean_Vol;
        else
            Ini_GBM_Vol   = -1;
            Ini_GARCH_Vol = -1;
        end
        Tenor = Tenors(ten);
        fprintf('\n Processing the %i-%s tenor...', Tenor, Settings.Time_Unit);    
        
        %% PAPER'S WAY OF SLICING DATA: Train #trajs windowed from start; Test #trajs windowed from back.
        % ----- slicing real DATA into n train trajs
        slicedTrainData = [];
        for t = (1 + Tenor) : T
            slicedTrainData    = [slicedTrainData; Underlying((t-Tenor):t-1)']; % moving window of 'Tennor' long, starting from 1. Results in 352x64 matrix (where 64 is tennor)
            [numTrainTrajs,cc] = size(slicedTrainData);
            if numTrainTrajs == Settings.numTrainTraj_realData(ten)
                break
            end
        end
        % figure,plot(slicedTrainData','x-'); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight;
        % title('Train Real Data, Sliced');
        % --- slicing real DATA into n test trajs
        slicedTestData = [];
        for t = T : -1 : 1
            slicedTestData    = [slicedTestData; Underlying((t-Tenor):t-1)']; % moving window of 'Tennor' long, starting from 1. Results in 352x64 matrix (where 64 is tennor)
            [numTestTrajs,cc] = size(slicedTestData);
            if numTestTrajs == Settings.numTestTraj_realData(ten)
                break
            end
        end
        % figure,plot(slicedTestData','x-'); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight; 
        % title('Test Real Data, Sliced')
        % ----- standardise all trajs to start from 1 i.e. divide all by S_1 (diff for individual trajs)
        % Train data
        [rr1,temp] = size(slicedTrainData);
        train_strikes = zeros(rr1,1);
        Data_Train = [];
        for ii = 1 : rr1
            train_strikes(ii,1) = slicedTrainData(ii,1); % train strikes
            Data_Train(ii,:)    = slicedTrainData(ii,:) / train_strikes(ii,1);
        end
        % figure,plot(Data_Train','x-'); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight;
        % title('PAPER Real Data Std Train')
        % Test data
        [rr2,temp] = size(slicedTestData);
        test_strikesRealData = zeros(rr2,1);
        Data_Test = [];
        for jj = 1 : rr2
            test_strikesRealData(jj,1) = slicedTestData(jj,1); % test strikes
            Data_Test(jj,:)            = slicedTestData(jj,:) / test_strikesRealData(jj,1);
        end
        % figure,plot(Data_Test','x-'); xlabel('Time [Days]'); ylabel('Stock Price [USD]'); axis tight;
        % title('PAPER Real Data Std Test')
        
        %% going over 3 RL algos: Settings.Methods = {'LSPI', 'FQI', 'LSM'};
        % --------------------- TRAIN ------------------------
        % initialWeights =  [0.8892; 0.0672; 0.7733; 0.2836; 0.2011; 0.6820; 0.7859];
        fprintf('\n Training phase \n');
        fprintf('\n Running the %s algorithm...', Settings.Methods{m});
        % ----- Train Weights with real data
        fprintf('\n Train weights with real data \n');
        initialWeights               = randn(numBasis,1);
        [W_Data,priceTrain_Data]     = Estimator(initialWeights, Data_Train, r,  Settings.Methods{m}, Realized_Vol_Window, GBM_Vol(ten,i), Settings);
        Weights{ten,i,1}             = W_Data;
        Prices_Train(ten,i,1)        = priceTrain_Data;
        % ----- Train Weights with GBM
        fprintf('\n Train weights with GBM \n');
        initialWeights               = randn(numBasis,1);
        [W_GBM,priceTrain_GBM]       = Estimator (initialWeights, GBM_Train,  r, Settings.Methods{m}, Realized_Vol_Window, GBM_Vol(ten,i), Settings);
        Weights{ten,i,2}             = W_GBM;
        Prices_Train(ten,i,2)        = priceTrain_GBM;
        % ----- Train Weights with GARCH
        fprintf('\n Train weights with GARCH \n');
        initialWeights               = randn(numBasis,1);
        [W_GARCH,priceTrain_GARCH]   = Estimator(initialWeights, GARCH_Train, r, Settings.Methods{m}, Realized_Vol_Window, GARCH_Mean_Vol, Settings);
        Weights{ten,i,3}             = W_GARCH;
        Prices_Train(ten,i,3)        = priceTrain_GARCH;
        
        % --------------- TEST DATA IS REAL DATA ---------------
        % --- Weights were trained with real data.
        method = 'Wdata, TestData';
        fprintf(['\n' method '\n'])
        [Ave_Payoff1,allPayoffs1,ExBoundary_root1,ExBoundary_exh1] = Performance(Settings.Methods{m}, W_Data,  Data_Test,   r, Realized_Vol_Window, Ini_GBM_Vol,   1, Settings, Tickers{i}, 0, method, K);
        Real_Payoffs(ten,i,1)        = Ave_Payoff1;
        Real_Payoffs_NonStd(ten,i,1) = Ave_Payoff1*K;
        Real_Boundaries{ten,i,1}     = ExBoundary_root1;
        Real_Boundaries_Exh{ten,i,1} = ExBoundary_exh1;
        % --- Weights were trained with data simulated with GBM.
        method = 'Wgbm, TestData';
        fprintf(['\n' method '\n'])
        [Ave_Payoff2,allPayoffs2,ExBoundary_root2,ExBoundary_exh2] = Performance(Settings.Methods{m}, W_GBM,   Data_Test,   r, Realized_Vol_Window, Ini_GBM_Vol,   1, Settings, Tickers{i}, 0, method, K);
        Real_Payoffs(ten,i,2)        = Ave_Payoff2;
        Real_Payoffs_NonStd(ten,i,2) = Ave_Payoff2*K;
        Real_Boundaries{ten,i,2}     = ExBoundary_root2;
        Real_Boundaries_Exh{ten,i,2} = ExBoundary_exh2;
        % --- Weights were trained with data simulated with GARCH.
        method = 'Wgarch, TestData';
        fprintf(['\n' method '\n'])
        [Ave_Payoff3,allPayoffs3,ExBoundary_root3,ExBoundary_exh3] = Performance(Settings.Methods{m}, W_GARCH, Data_Test,   r, Realized_Vol_Window, Ini_GARCH_Vol, 1, Settings, Tickers{i}, 0, method, K);
        Real_Payoffs(ten,i,3)        = Ave_Payoff3;
        Real_Payoffs_NonStd(ten,i,3) = Ave_Payoff3*K;
        Real_Boundaries{ten,i,3}     = ExBoundary_root3;
        Real_Boundaries_Exh{ten,i,3} = ExBoundary_exh3;
      
        % --------------- TEST DATA IS GBM AND GARCH SIMULATED DATA ---------------
        % --- Weights were trained with GBM simulated data.
        method = 'Wgbm, TestGbm';
        fprintf(['\n' method '\n'])
        [Ave_Payoff4,allPayoffs4,ExBoundary_root4,ExBoundary_exh4] = Performance(Settings.Methods{m}, W_GBM,    GBM_Test,   r, Realized_Vol_Window, Ini_GBM_Vol,   1, Settings, Tickers{i}, 1, method, K);
        Sim_Payoffs(ten,i,1)         = Ave_Payoff4;
        Sim_Payoffs_NonStd(ten,i,1)  = Ave_Payoff4*K;
        Sim_Boundaries{ten,i,1}      = ExBoundary_root4;
        Sim_Boundaries_exh{ten,i,1}  = ExBoundary_exh4;
        % --- Weights were trained with GARCH simulated data.
        method = 'Wgarch, TestGarch';
        fprintf(['\n' method '\n'])
        [Ave_Payoff5,allPayoffs5,ExBoundary_root5,ExBoundary_exh5]  = Performance(Settings.Methods{m}, W_GARCH, GARCH_Test, r, Realized_Vol_Window, Ini_GARCH_Vol, 1, Settings, Tickers{i}, 1, method, K);
        Sim_Payoffs(ten,i,2)         = Ave_Payoff5;
        Sim_Payoffs_NonStd(ten,i,2)  = Ave_Payoff5*K;
        Sim_Boundaries{ten,i,2}      = ExBoundary_root5;
        Sim_Boundaries_exh{ten,i,2}  = ExBoundary_exh5;
        fprintf(['\n The %i-%s tenor has been processed.', Tenor, Settings.Time_Unit, '\n']);
    end % end tenor
end % end Stocks I

fprintf(['The %s algorithm has been run.', Settings.Methods{m}, '\n']);

OutputName                  = Settings.Methods{m}; % Algo name
results.GBM_Vol             = GBM_Vol;
results.Weights             = Weights;
results.Prices_Train        = Prices_Train;
results.Real_Payoffs        = Real_Payoffs;
results.Real_Boundaries     = Real_Boundaries;
results.Real_Boundaries_Exh = Real_Boundaries_Exh;
results.Real_Payoffs_NonStd = Real_Payoffs_NonStd;
results.Sim_Payoffs         = Sim_Payoffs;
results.Sim_Payoffs_NonStd  = Sim_Payoffs_NonStd;
results.Sim_Boundaries      = Sim_Boundaries;
results.Sim_Boundaries_exh  = Sim_Boundaries_exh;
save(['results/' OutputName '.mat'], 'results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAYING THE RESULTS





fprintf( '\n\n\nAVERAGE PAYOFFS ON REAL DATA:\n' );
% % fprintf( ['\n',FWS('',12,'left')] );
% % for m = 1:3
% %     fprintf( [RPS(' ',10), FWS(Settings.Methods{m},10,'right'), RPS(' ',10)] );
% % end
% % fprintf( ['\n',FWS('Maturity',12,'left')] );
% % for m = 1:3
% %     fprintf( [FWS('GBM',10,'right'), FWS('GARCH',10,'right'), FWS('Data',10,'right')] );
% % end
% % fprintf( ['\n',RPS('-',12 + 10 * 3 * 3)] );
% % for ten = 1:TEN
% %     fprintf( ['\n', FWS(sprintf('%i %ss',Tenors(ten),Settings.Time_Unit),12,'left')] );
% %
% %     for m = 1:3
% %         for ds = 1:3
% %             ind         = isnan(Real_Payoffs(ten,m,ds,:));
% %             fprintf( FWS(sprintf('%2.3f',mean(Real_Payoffs(ten,m,ds,~ind))),10,'right') );
% %         end
% %     end
% % end
% % fprintf( ['\n',RPS('-',12 + 10 * 3 * 3)] );
% %
% % fprintf( '\n\n\nAVERAGE PAYOFFS ON SIMULATED DATA:\n' );
% % fprintf( ['\n',FWS('',12,'left')] );
% % for m = 1:3
% %     fprintf( [RPS(' ',11), FWS(Settings.Methods{m},4,'right'), RPS(' ',5)] );
% % end
% % fprintf( ['\n',FWS('Maturity',12,'left')] );
% % for m = 1:3
% %     fprintf( [FWS('GBM',10,'right'), FWS('GARCH',10,'right')] );
% % end
% % fprintf( ['\n',RPS('-',12 + 10 * 3 * 2)] );
% % for ten = 1:TEN
% %     fprintf( ['\n', FWS(sprintf('%i %ss',Tenors(ten),Settings.Time_Unit),12,'left')] );
% %
% %     for m = 1:3
% %         for ds = 1:2
% %             ind         = isnan(Sim_Payoffs(ten,m,ds,:));
% %             fprintf( FWS(sprintf('%2.3f',mean(Sim_Payoffs(ten,m,ds,~ind))),10,'right') );
% %         end
% %     end
% % end
% % fprintf( ['\n',RPS('-',12 + 10 * 3 * 2)] );
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % T-TESTS
% %
% % fprintf( '\n\n\nP-VALUES OF PAIRED SAMPLE T-TESTS FOR AVERAGE PAYOFFS:\n' );
% % for ten = 1:TEN
% %     fprintf('\n   For the %i-%s tenor:', Tenors(ten), Settings.Time_Unit);
% %
% %     for m1 = 1:3
% %         for m2 = 1:3
% %             if m1 < m2
% %                 fprintf('\n      %s vs %s:', Settings.Methods{m1}, Settings.Methods{m2});
% %                 [h,p_GBM]       = ttest(Real_Payoffs(ten,m1,1,:), Real_Payoffs(ten,m2,1,:));
% %                 [h,p_GARCH]     = ttest(Real_Payoffs(ten,m1,2,:), Real_Payoffs(ten,m2,2,:));
% %                 [h,p_Data]      = ttest(Real_Payoffs(ten,m1,3,:), Real_Payoffs(ten,m2,3,:));
% %                 fprintf('\n         Real / simulated training data & real test data: GBM = %0.5f, GARCH = %0.5f, real data = %0.5f.', p_GBM, p_GARCH, p_Data);
% %                 [h,p_GBM]       = ttest(Sim_Payoffs(ten,m1,1,:), Sim_Payoffs(ten,m2,1,:));
% %                 [h,p_GARCH]     = ttest(Sim_Payoffs(ten,m1,2,:), Sim_Payoffs(ten,m2,2,:));
% %                 fprintf('\n         Simulated training data & simulated test data: GBM = %0.5f, GARCH = %0.5f.', p_GBM, p_GARCH);
% %             end
% %         end
% %     end
% %     fprintf('\n');
% %     for m = 1:3
% %         fprintf('\n      For %s:', Settings.Methods{m});
% %         [h,p_GBM_GARCH]     = ttest(Real_Payoffs(ten,m,1,:), Real_Payoffs(ten,m,2,:));
% %         [h,p_GBM_Data]      = ttest(Real_Payoffs(ten,m,1,:), Real_Payoffs(ten,m,3,:));
% %         [h,p_GARCH_Data]    = ttest(Real_Payoffs(ten,m,2,:), Real_Payoffs(ten,m,3,:));
% %         fprintf('\n         Real / simulated training data & real test data: GBM vs GARCH = %0.5f, GBM vs real data = %0.5f, GARCH vs real data = %0.5f.', p_GBM_GARCH, p_GBM_Data, p_GARCH_Data);
% %         [h,p_GBM_GARCH]     = ttest(Sim_Payoffs(ten,m,1,:), Sim_Payoffs(ten,m,2,:));
% %         fprintf('\n         Simulated training data & simulated test data: GBM vs GARCH = %0.5f.', p_GBM_GARCH);
% %     end
% %     fprintf('\n\n');
% % end