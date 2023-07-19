classdef Xavi_Box < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        VisualizationSwitch             matlab.ui.control.Switch
        VisualizationSwitchLabel        matlab.ui.control.Label
        zoffsetinmmEditField            matlab.ui.control.NumericEditField
        OffsetinmmLabel                 matlab.ui.control.Label
        UnitButtonGroup                 matlab.ui.container.ButtonGroup
        mmButton                        matlab.ui.control.RadioButton
        mButton                         matlab.ui.control.RadioButton
        inchButton                      matlab.ui.control.RadioButton
        WirelengthEditField             matlab.ui.control.EditField
        WirelengthEditFieldLabel        matlab.ui.control.Label
        InductanceEditField             matlab.ui.control.EditField
        InductanceEditFieldLabel        matlab.ui.control.Label
        WorstcoilwireresolutionEditField  matlab.ui.control.EditField
        WorstcoilwireresolutionEditFieldLabel  matlab.ui.control.Label
        AvgcoilwireresolutionEditField  matlab.ui.control.EditField
        AvgcoilwireresolutionEditFieldLabel  matlab.ui.control.Label
        CoordinateSystemSwitch          matlab.ui.control.Switch
        CoordinateSystemLabel           matlab.ui.control.Label
        KodeinBoxTMSKCoilDesignInstrumentLabel  matlab.ui.control.Label
        SaveNIfTIButton                 matlab.ui.control.Button
        CalculateAFieldButton           matlab.ui.control.Button
        LoadCoilConductorPathButton     matlab.ui.control.Button
        UIAxes                          matlab.ui.control.UIAxes
    end

%%
%   TMS KoDeIn Box, Version 1.0, 12/2022
%   Kaiserslautern Coil Design and Instrument Box
%   Max Koehler and Stefan M. Goetz
%   (c) 2006-2022



    properties (Access = public)
        coil
        coilraw
        dAdt
        niftistring
        dAdtcoil
        vol
        sup
        xx
        yy
        zz
        AFlag
        hdr
    end
    
    methods (Access = private)
%%  Function Update User Interface

        function updateui(app)

            for n=1:size(app.coilraw,3)

                plot3(app.UIAxes,app.coil(:,1,n),app.coil(:,2,n),app.coil(:,3,n),'k')
                
                if n==1
                    hold (app.UIAxes,'on')
                end

            end
            
            hold (app.UIAxes,'off')
            
            if app.AFlag ==1

                if app.VisualizationSwitch.Value == 0

                    hold (app.UIAxes,'on')
                    
                    %// Get the current colormap
                    currentColormap = colormap(hot);
                    close(gcf)
                    
                    n=5; %set resolution
                
                    q=quiver3(app.UIAxes,app.xx(1:n:end,1:n:end,1:n:end),app.yy(1:n:end,1:n:end,1:n:end),app.zz(1:n:end,1:n:end,1:n:end),app.vol(1:n:end,1:n:end,1:n:end,2),app.vol(1:n:end,1:n:end,1:n:end,1),app.vol(1:n:end,1:n:end,1:n:end,3),2.5);
                
                    %// Compute the magnitude of the vectors
                    mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
                            reshape(q.WData, numel(q.UData), [])).^2, 2));
                
                    %// Now determine the color to make each arrow using a colormap
                    [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
                
                    %// Now map this to a colormap to get RGB
                    cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
                    cmap(:,:,4) = 255;
                    cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
                
                    %// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
                    set(q.Head, ...
                        'ColorBinding', 'interpolated', ...
                        'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
                
                    %// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
                    set(q.Tail, ...
                        'ColorBinding', 'interpolated', ...
                        'ColorData', reshape(cmap(1:2,:,:), [], 4).');
    
                    hold (app.UIAxes,'off')
                     
                elseif app.VisualizationSwitch.Value == 1

                    hold (app.UIAxes,'on')
                    n=1;
                    mysplstd = 100;
                    spts = 400;
                    
                    ccenterx = mean(app.coil(:,1));
                    ccentery = mean(app.coil(:,2));
                    ccenterz = mean(app.coil(:,3));
                    startxx = ccenterx + mysplstd* randn(1, spts);
                    startyy = ccentery + mysplstd* randn(1, spts); 
                    startzz = ccenterz + mysplstd* randn(1, spts);
                    
                    mysl = streamline(app.UIAxes,stream3(app.xx(1:n:end,1:n:end,1:n:end), app.yy(1:n:end,1:n:end,1:n:end), app.zz(1:n:end,1:n:end,1:n:end), app.vol(1:n:end,1:n:end,1:n:end,2), app.vol(1:n:end,1:n:end,1:n:end,1), app.vol(1:n:end,1:n:end,1:n:end,3), startxx, startyy, startzz));
                    hold (app.UIAxes,'off')
                    
                end
            end
            

            for n=1:size(app.coilraw,3)

                cdiff = diff(app.coil(:,:,n));
                ncdiff(:,n) = sqrt(cdiff(:,1).^2+cdiff(:,2).^2+cdiff(:,3).^2);

            end

            maxdiff = max(max(ncdiff));
            mncdiff = mean(reshape(ncdiff,[],1));

            b1=6.5;
            b2=20;

            if mncdiff < b1
                app.AvgcoilwireresolutionEditField.BackgroundColor = "green";
            elseif (b1 <= mncdiff) && (mncdiff <= b2)
                app.AvgcoilwireresolutionEditField.BackgroundColor = "yellow";
            elseif b2 < mncdiff
                app.AvgcoilwireresolutionEditField.BackgroundColor = "red";
            end

            if maxdiff < b1
                app.WorstcoilwireresolutionEditField.BackgroundColor = "green";
            elseif (b1 <= maxdiff) && (maxdiff <= b2)
                app.WorstcoilwireresolutionEditField.BackgroundColor = "yellow";
            elseif b2 < maxdiff
                app.WorstcoilwireresolutionEditField.BackgroundColor = "red";
            end

            app.AvgcoilwireresolutionEditField.Value=[num2str(mncdiff,'%.2f'),' mm'];
            app.WorstcoilwireresolutionEditField.Value=[num2str(maxdiff,'%.2f'),' mm'];

        end
        
%% Function Calculate Inductance and Wire Length        
        function induct_wireleng(app)
    
            wireradius = 0.000125;                        % Radius of the wire in meters (on axis towards origin) Default is in mm = 2.5e-3   
            cdistribfactor = 1/3;                       % Current distribution (0: perfect skin effect, 1/4: equal/homogeneous distrib)
            mu0 = 4*pi*1E-7;
            
            for n=1:size(app.coilraw,3)

                coilm=app.coil(:,:,n).*1e-5;% -default = -3    3 when the wire radius is in (1cm), -4 for mm radius, -5 for micrometer radius
    
                Pdiff = diff(coilm);
                mylength(n) = nansum(sqrt(sum(Pdiff.^2, 2)));
    
                section=mylength(n)/size(coilm,1);
                
                addpoints=round(norm(coilm(1,:)-coilm(end,:))/section);
    
                pts(:,1)=linspace(coilm(end,1),coilm(1,1),addpoints+1);
                pts(:,2)=linspace(coilm(end,2),coilm(1,2),addpoints+1);
                pts(:,3)=linspace(coilm(end,3),coilm(1,3),addpoints+1);
    
                pts=pts(2:end-1,:);
    
                coilm=[coilm;pts];
                
                coilseg(:,1) = interp1(coilm(:,1), 1:1/2:size(coilm,1), 'pchip');
                coilseg(:,2) = interp1(coilm(:,2), 1:1/2:size(coilm,1), 'pchip');
                coilseg(:,3) = interp1(coilm(:,3), 1:1/2:size(coilm,1), 'pchip');
                
                % init:
                Lfact = 0;
                
                for icnt=2:size(coilseg,1)
                        
                    D = sqrt((coilseg(icnt,1)-coilseg(2:end,1)).^2 + (coilseg(icnt,2)-coilseg(2:end,2)).^2 + (coilseg(icnt,3)-coilseg(2:end,3)).^2);    % starting with 2:end, because need for Lfact the same size (N points lead to N-1 segments!)
                    Pdiff = diff(coilseg);
                    Lfact = Lfact + sum((((coilseg(icnt,1)-coilseg(icnt-1,1))*Pdiff(:,1)) + ((coilseg(icnt,2)-coilseg(icnt-1,2))*Pdiff(:,2)) + ((coilseg(icnt,3)-coilseg(icnt-1,3))*Pdiff(:,3)))./max(D, wireradius));
                    
                end
                
                L(n) = mu0/(4*pi) * (Lfact + 2*mylength(n)*cdistribfactor);
                
                clear section addpoints pts Lfact coilm coilseg

            end

            av_L=mean(L);
            av_mylength=mean(mylength);

            if av_L > 0
                app.InductanceEditField.BackgroundColor = "green";
            else
                app.InductanceEditField.BackgroundColor = "red";
            end

            if av_mylength ~= 0
                app.WirelengthEditField.BackgroundColor = "green";
            else
                app.WirelengthEditField.BackgroundColor = "red";
            end

            app.InductanceEditField.Value=[num2str(av_L*1e6,'%.2f'),' µH']; %defalut (av_L*1e6,'%.2f')
            app.WirelengthEditField.Value=[num2str(av_mylength,'%.2f'),' m'];
        end
    
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadCoilConductorPathButton
        function LoadCoilConductorPathButtonPushed(app, event)
            
            app.AFlag=0;

            [myfilename, mypath]=uigetfile({'*.mat; *.xlsx; *.xls; *.csv','Supported coil wire path files (*.mat, *.xlsx, *.xls, *.csv)'; '*.mat','Matlab data files (*.mat)'; '*.xlsx; *.xls','Excel file (*.xlsx, *.xls)'; '*.csv','Comma-separated values (*.csv)'}, 'Select coil conductor path file ...');
            
            if (~isequal(myfilename, 0)) && (~isequal(mypath, 0))    % both 0 means user clicked on cancel
            
                % check format
                if strfind(myfilename(length(myfilename)-3:length(myfilename)), '.mat')
                    % Matlab .mat file
                
                    app.coilraw=importdata([mypath myfilename]);          
                
                elseif strfind(myfilename(length(myfilename)-4:length(myfilename)), '.xlsx')
                    % xml Excel file
                    % read from line two on
                    
                    tmp_xlsfilesheet1 = xlsread([mypath myfilename]);
                    app.coilraw = tmp_xlsfilesheet1;
                    
                    clear tmp_xlsfilesheet1
                    
                elseif strfind(myfilename(length(myfilename)-3:length(myfilename)), '.xls')
                    % old Excel file
                    % read from line two on
                    
                    tmp_xlsfilesheet1 = xlsread([mypath myfilename]);
                    app.coilraw = tmp_xlsfilesheet1; 
                    
                    clear tmp_xlsfilesheet1
                    
                elseif strfind(myfilename(length(myfilename)-3:length(myfilename)), '.csv')
                    % CSV file
                    % read from lie two on
                    
                    tmp_csvfile = importdata([mypath myfilename]);
                    app.coilraw = tmp_csvfile.data;
                    
                    clear tmp_csvfile
                    
                else
                    % unknown file type
                    msgbox('Unknown file type.', 'Coil file', 'error')
                    return
                end
                
                clear oldFolder myfilename mypath

                for  n=1:fix(size(app.coilraw,2)/3)

                    coilrawtmp(:,:,n)=app.coilraw(:,1+((n-1)*3):3+((n-1)*3));

                end

                app.coilraw=coilrawtmp;

                app.coil=app.coilraw;
                
                UnitButtonGroupSelectionChanged(app)
                induct_wireleng(app);
                updateui(app);
                
            else
                % in case they clicked on cancel
                
            end
                
        end

        % Selection changed function: UnitButtonGroup
        function UnitButtonGroupSelectionChanged(app, event)
            
            selectedButton = app.UnitButtonGroup.SelectedObject;
            value = app.zoffsetinmmEditField.Value;
            
            if app.AFlag == 0

                if app.mmButton.Value == true
                    app.coil=app.coilraw;
                elseif app.mButton.Value == true
                    app.coil=app.coilraw.*1000;
                elseif app.inchButton.Value == true
                    app.coil=app.coilraw.*25.4;
                end
    
                if app.mmButton.Value == true
                    app.coil(:,3,:)=app.coilraw(:,3,:)+value;
                elseif app.mButton.Value == true
                    app.coil(:,3,:)=(app.coilraw(:,3,:)+(value/1000)).*1000;
                elseif app.inchButton.Value == true
                    app.coil(:,3,:)=(app.coilraw(:,3,:)+(value/25.4)).*25.4;
                end
    
                induct_wireleng(app);
                updateui(app);
    
                if app.inchButton.Value == true
                    [y,Fs] = audioread('eagle.mp3');
                    player=audioplayer(y,Fs);
                    playblocking(player)
                end

            end

        end

        % Value changed function: zoffsetinmmEditField
        function zoffsetinmmEditFieldValueChanged(app, event)
            
            value = app.zoffsetinmmEditField.Value;
            
            if app.mmButton.Value == true
                app.coil(:,3,:)=app.coilraw(:,3,:)+value;
            elseif app.mButton.Value == true
                app.coil(:,3,:)=(app.coilraw(:,3,:)+(value/1000)).*1000;
            elseif app.inchButton.Value == true
                app.coil(:,3,:)=(app.coilraw(:,3,:)+(value/25.4)).*25.4;
            end

            induct_wireleng(app);
            updateui(app);

        end

        % Value changed function: CoordinateSystemSwitch
        function CoordinateSystemSwitchValueChanged(app, event)
            
            if app.AFlag == 0

                app.coil(:,[1 2 3],:)=app.coil(:,[2 1 3],:);
                updateui(app);

            end

        end

        % Value changed function: VisualizationSwitch
        function VisualizationSwitchValueChanged(app, event)
            
            value = app.VisualizationSwitch.Value;
            
            if app.AFlag == 1

                updateui(app);

            end

        end

        % Button pushed function: CalculateAFieldButton
        function CalculateAFieldButtonPushed(app, event)
            
            addpath('simnibs_functions');
    	    app.hdr = nifti_load('nifti_blank.nii');
            
            app.hdr.sform(3,4)=-85;

            [app.xx,app.yy,app.zz]=meshgrid(linspace(-((app.hdr.dim(2)-1)*app.hdr.sform(1,1)+app.hdr.sform(1,4)),(app.hdr.dim(2)-1)*app.hdr.sform(1,1)+app.hdr.sform(1,4),app.hdr.dim(2)),...
                                            linspace(-((app.hdr.dim(3)-1)*app.hdr.sform(2,2)+app.hdr.sform(2,4)),(app.hdr.dim(3)-1)*app.hdr.sform(2,2)+app.hdr.sform(2,4),app.hdr.dim(3)),...
                                            linspace((app.hdr.dim(4)-1)*app.hdr.sform(3,3)+app.hdr.sform(3,4),app.hdr.sform(3,4),app.hdr.dim(4)));

            app.sup(:,1)=reshape(app.xx,[],1);
            app.sup(:,2)=reshape(app.yy,[],1);
            app.sup(:,3)=reshape(app.zz,[],1);

            app.vol=zeros(size(app.xx,1),size(app.xx,2),size(app.xx,3),3)

            app.hdr.sform(3,4)=0;

            Mu0=            1E-7;           % Magnetic field constant/4*pi
            preconst=       Mu0;

            sups=size(app.sup,1);
            supi=app.sup;

            d = uiprogressdlg(app.UIFigure,'Icon','membrane.png', 'Title','Please Wait','Message','Calculating A-Field');
            d.Value = 0;

            for n=1:size(app.coilraw,3)
                
                coil2=app.coil(:,:,n);
                coilvec=diff(coil2);
                coil2(end,:)=[];

                for i=1:sups    % parfor in original
                    
                    normdist = sqrt((coil2(:,1)-supi(i,1)).*(coil2(:,1)-supi(i,1))...
                                   +(coil2(:,2)-supi(i,2)).*(coil2(:,2)-supi(i,2))...
                                   +(coil2(:,3)-supi(i,3)).*(coil2(:,3)-supi(i,3)));
            
                    dAdtcoil(i,:) = preconst*sum((1./normdist).*coilvec(:,:),1);
            
                end
                
                app.vol(:,:,:,2)=app.vol(:,:,:,2)+reshape(dAdtcoil(:,1),size(app.yy,1),size(app.yy,2),size(app.yy,3))/size(app.coilraw,3); %divide by the number of strands
                app.vol(:,:,:,1)=app.vol(:,:,:,1)+reshape(dAdtcoil(:,2),size(app.xx,1),size(app.xx,2),size(app.xx,3))/size(app.coilraw,3);
                app.vol(:,:,:,3)=app.vol(:,:,:,3)+reshape(dAdtcoil(:,3),size(app.zz,1),size(app.zz,2),size(app.zz,3))/size(app.coilraw,3);

                clear dAdtcoil   

                d.Value = d.Value + 1/size(app.coilraw,3);

            end

            app.AFlag=1;
            close(d)
            updateui(app);
            app.hdr.vol=app.vol;

        end

        % Button pushed function: SaveNIfTIButton
        function SaveNIfTIButtonPushed(app, event)
            
            [file, path] = uiputfile({'*.nii','NIfTI File (*.nii)'}); 
            
            oldFolder = cd;
            cd (path);
            nifti_save(app.hdr,file);
            cd (oldFolder);

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 745 500];
            app.UIFigure.Name = 'Coil Path NIfTI Generator';
            app.UIFigure.Icon = fullfile(pathToMLAPP, 'LitzIcon_transparent.png');

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            xlabel(app.UIAxes, 'x (mm)')
            ylabel(app.UIAxes, 'y (mm)')
            zlabel(app.UIAxes, 'z (mm)')
            app.UIAxes.View = [22.5 45];
            app.UIAxes.DataAspectRatio = [1 1 1];
            app.UIAxes.PlotBoxAspectRatio = [1.35593220338983 1.35593220338983 1];
            app.UIAxes.XLim = [-200 200];
            app.UIAxes.YLim = [-200 200];
            app.UIAxes.ZLim = [-85 210];
            app.UIAxes.XAxisLocation = 'origin';
            app.UIAxes.XTick = [-200 0 200];
            app.UIAxes.XTickLabel = {'-200'; '0'; '200'};
            app.UIAxes.YAxisLocation = 'origin';
            app.UIAxes.YTick = [-200 0 200];
            app.UIAxes.YTickLabel = {'-200'; '0'; '200'};
            app.UIAxes.ZTick = [-85 0 210];
            app.UIAxes.ZTickLabel = {'-85'; '0'; '210'};
            app.UIAxes.BoxStyle = 'full';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.ZGrid = 'on';
            app.UIAxes.ColorOrder = [0 0.451 0.7412;0.851 0.3255 0.098;0.9294 0.6941 0.1255;0.4941 0.1843 0.5569;0.4667 0.6745 0.1882;0.302 0.7451 0.9333;0.6353 0.0784 0.1843];
            app.UIAxes.Position = [296 68 434 374];

            % Create LoadCoilConductorPathButton
            app.LoadCoilConductorPathButton = uibutton(app.UIFigure, 'push');
            app.LoadCoilConductorPathButton.ButtonPushedFcn = createCallbackFcn(app, @LoadCoilConductorPathButtonPushed, true);
            app.LoadCoilConductorPathButton.Tooltip = {'Load file (mat, excel, comma-separated-value file) of a coil path or several parallel paths in millimeters, meters, or inchs (set unit below after loading the file). If the coil is not correctly displayed in the graph on the right after loading the file, adjust the units and/or the z offset.'};
            app.LoadCoilConductorPathButton.Position = [78 405 153 37];
            app.LoadCoilConductorPathButton.Text = 'Load Coil Conductor Path';

            % Create CalculateAFieldButton
            app.CalculateAFieldButton = uibutton(app.UIFigure, 'push');
            app.CalculateAFieldButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateAFieldButtonPushed, true);
            app.CalculateAFieldButton.Tooltip = {'Calculate the magnetic vector potential A for the coil path. Calculations might take a while. The results will be displayed on the right.'};
            app.CalculateAFieldButton.Position = [101 52 106 22];
            app.CalculateAFieldButton.Text = 'Calculate A-Field';

            % Create SaveNIfTIButton
            app.SaveNIfTIButton = uibutton(app.UIFigure, 'push');
            app.SaveNIfTIButton.ButtonPushedFcn = createCallbackFcn(app, @SaveNIfTIButtonPushed, true);
            app.SaveNIfTIButton.Tooltip = {'Save the calculated magnetic vector potential to a file of your choosing to be used in field modeling tools afterwards.'};
            app.SaveNIfTIButton.Position = [104 17 100 22];
            app.SaveNIfTIButton.Text = 'Save NIfTI';

            % Create KodeinBoxTMSKCoilDesignInstrumentLabel
            app.KodeinBoxTMSKCoilDesignInstrumentLabel = uilabel(app.UIFigure);
            app.KodeinBoxTMSKCoilDesignInstrumentLabel.FontSize = 24;
            app.KodeinBoxTMSKCoilDesignInstrumentLabel.Position = [36 456 545 32];
            app.KodeinBoxTMSKCoilDesignInstrumentLabel.Text = 'Kodein Box – TMS (K)Coil Design Instrument';

            % Create CoordinateSystemLabel
            app.CoordinateSystemLabel = uilabel(app.UIFigure);
            app.CoordinateSystemLabel.HorizontalAlignment = 'center';
            app.CoordinateSystemLabel.Position = [107 224 107 28];
            app.CoordinateSystemLabel.Text = 'Coordinate System';

            % Create CoordinateSystemSwitch
            app.CoordinateSystemSwitch = uiswitch(app.UIFigure, 'slider');
            app.CoordinateSystemSwitch.Items = {'left', 'right'};
            app.CoordinateSystemSwitch.ValueChangedFcn = createCallbackFcn(app, @CoordinateSystemSwitchValueChanged, true);
            app.CoordinateSystemSwitch.Tooltip = {'Switch between right-hand and left-hand (often used in Matlab) coordinate systems, essentially swapping x and y.'};
            app.CoordinateSystemSwitch.Position = [123 258 54 24];
            app.CoordinateSystemSwitch.Value = 'left';

            % Create AvgcoilwireresolutionEditFieldLabel
            app.AvgcoilwireresolutionEditFieldLabel = uilabel(app.UIFigure);
            app.AvgcoilwireresolutionEditFieldLabel.Position = [46 193 131 22];
            app.AvgcoilwireresolutionEditFieldLabel.Text = 'Avg. coil wire resolution';

            % Create AvgcoilwireresolutionEditField
            app.AvgcoilwireresolutionEditField = uieditfield(app.UIFigure, 'text');
            app.AvgcoilwireresolutionEditField.HorizontalAlignment = 'center';
            app.AvgcoilwireresolutionEditField.Tooltip = {'Spacing of sample points of the coil path representing the resolution of the coil description.'};
            app.AvgcoilwireresolutionEditField.Position = [185 193 77 22];

            % Create WorstcoilwireresolutionEditFieldLabel
            app.WorstcoilwireresolutionEditFieldLabel = uilabel(app.UIFigure);
            app.WorstcoilwireresolutionEditFieldLabel.Position = [46 157 138 22];
            app.WorstcoilwireresolutionEditFieldLabel.Text = 'Worst coil wire resolution';

            % Create WorstcoilwireresolutionEditField
            app.WorstcoilwireresolutionEditField = uieditfield(app.UIFigure, 'text');
            app.WorstcoilwireresolutionEditField.HorizontalAlignment = 'center';
            app.WorstcoilwireresolutionEditField.Position = [185 157 77 22];

            % Create InductanceEditFieldLabel
            app.InductanceEditFieldLabel = uilabel(app.UIFigure);
            app.InductanceEditFieldLabel.Position = [46 122 96 22];
            app.InductanceEditFieldLabel.Text = 'Inductance';

            % Create InductanceEditField
            app.InductanceEditField = uieditfield(app.UIFigure, 'text');
            app.InductanceEditField.HorizontalAlignment = 'center';
            app.InductanceEditField.Tooltip = {'Magnetic inductance of the coil (without any cable or connectors attached). Most TMS devices expect around 10 µH.'};
            app.InductanceEditField.Position = [185 122 77 22];

            % Create WirelengthEditFieldLabel
            app.WirelengthEditFieldLabel = uilabel(app.UIFigure);
            app.WirelengthEditFieldLabel.Position = [46 87 96 22];
            app.WirelengthEditFieldLabel.Text = 'Wire length';

            % Create WirelengthEditField
            app.WirelengthEditField = uieditfield(app.UIFigure, 'text');
            app.WirelengthEditField.HorizontalAlignment = 'center';
            app.WirelengthEditField.Tooltip = {'Required length of the conductor for the coil only if the coil is to be implemented (ignoring the coil cable).'};
            app.WirelengthEditField.Position = [185 87 77 22];

            % Create UnitButtonGroup
            app.UnitButtonGroup = uibuttongroup(app.UIFigure);
            app.UnitButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @UnitButtonGroupSelectionChanged, true);
            app.UnitButtonGroup.Tooltip = {'You may have to switch between the units to adjust them to your coil path file. It often is the solution if you do not see a coil in the graph on the side after loading a file.'};
            app.UnitButtonGroup.TitlePosition = 'centertop';
            app.UnitButtonGroup.Title = 'Unit';
            app.UnitButtonGroup.Position = [59 336 190 55];

            % Create inchButton
            app.inchButton = uiradiobutton(app.UnitButtonGroup);
            app.inchButton.Text = 'inch';
            app.inchButton.Position = [134 7 48 22];

            % Create mButton
            app.mButton = uiradiobutton(app.UnitButtonGroup);
            app.mButton.Text = 'm';
            app.mButton.Position = [72 7 48 22];

            % Create mmButton
            app.mmButton = uiradiobutton(app.UnitButtonGroup);
            app.mmButton.Text = 'mm';
            app.mmButton.Position = [9 7 48 22];
            app.mmButton.Value = true;

            % Create OffsetinmmLabel
            app.OffsetinmmLabel = uilabel(app.UIFigure);
            app.OffsetinmmLabel.HorizontalAlignment = 'center';
            app.OffsetinmmLabel.Position = [82 300 87 22];
            app.OffsetinmmLabel.Text = 'z offset (in mm)';

            % Create zoffsetinmmEditField
            app.zoffsetinmmEditField = uieditfield(app.UIFigure, 'numeric');
            app.zoffsetinmmEditField.Limits = [-200 200];
            app.zoffsetinmmEditField.ValueChangedFcn = createCallbackFcn(app, @zoffsetinmmEditFieldValueChanged, true);
            app.zoffsetinmmEditField.Tooltip = {'Adjust the zero offset in z direction, e.g., to 115 mm for the zero position of SimNIBS or to add the thickness of the bottom shell of a coiil and account for the thickness of the coil conductor.'};
            app.zoffsetinmmEditField.Position = [171 300 56 22];

            % Create VisualizationSwitchLabel
            app.VisualizationSwitchLabel = uilabel(app.UIFigure);
            app.VisualizationSwitchLabel.HorizontalAlignment = 'center';
            app.VisualizationSwitchLabel.Position = [476 3 72 22];
            app.VisualizationSwitchLabel.Text = 'Visualization';

            % Create VisualizationSwitch
            app.VisualizationSwitch = uiswitch(app.UIFigure, 'slider');
            app.VisualizationSwitch.Items = {'Vector', 'Stream'};
            app.VisualizationSwitch.ItemsData = [0 1];
            app.VisualizationSwitch.ValueChangedFcn = createCallbackFcn(app, @VisualizationSwitchValueChanged, true);
            app.VisualizationSwitch.Tooltip = {'Switch between vector arrows and stream line illustration of the magnetic vector potential A.'};
            app.VisualizationSwitch.Position = [488 34 45 20];
            app.VisualizationSwitch.Value = 0;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Xavi_Box

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end