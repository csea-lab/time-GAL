function circularConnectivity125times(MatrixData, lengthTime, timestep, legends, varargin)
%   Draws the circular graph for 125 channels from EGI EEG system
%   The function is fed with the GAL results and the mask of statistically
%   significant points. This function is used  for each time step of the
%   trial.
%
%   Input:
%       - MatrixData : Matrix of GAL data with all decoding rate results
%       averaged across all subjects.
%       - lengthTime : Length of the time step, also known as delta Time.
%       - timestep : Indicates the moment or point of time off this plot.
%       - legends : 0 or 1 to show the bottom legend of what dots
%       represent.
%   
%   MatrixData must be square matrices with same number of rows and
%   columns.

for var = 1:length(varargin)
    if strcmp(varargin(var),'Gradient') % switch if more args
        gradient = true;
    else
        gradient = false;
    end
end


%%  1. Initialize values and draw EEG circunference

    hold on
    
    % Legend of Channels =  1 FL; 2 CL; 3 TL; 4 PL; 5 OL; 6 OR; 7 PR; 8 TR; 9 CR; 10 FR
    infoCh = [10 10 10 10 10 2 2 10 10 10 10 1 2 10 10 1 1 1 1 1 1 1 ...
        1 1 1 1 1 1 2 2 2 1 1 1 2 2 2 3 3 3 2 2 3 3 3 3 2 3 3 3 3 4 4 4 4 3 ...
        3 3 4 4 4 4 3 3 4 4 4 5 5 5 5 7 5 5 5 6 7 7 7 9 6 6 6 7 7 7 9 6 6 7 ...
        7 7 9 6 8 8 8 9 8 8 8 8 9 9 9 9 8 8 8 9 9 9 8 8 8 10 10 10 8 8 8 10 10 10 7];
    % reorder by groups
    FL = find(infoCh == 1); CL = find(infoCh == 2); TL = find(infoCh == 3); PL = find(infoCh == 4); OL = find(infoCh == 5);
    FR = find(infoCh == 10); CR = find(infoCh == 9); TR = find(infoCh == 8); PR = find(infoCh == 7); OR = find(infoCh == 6);
    
    %Colors of the channels groups
    Colors = [142, 117, 255; 255, 126, 117; 204, 255, 117; 255,164,58; 117,158,255;  117,158,255; 255,164,58; 227, 255, 117; 255, 126, 117; 142, 117, 255]./255;
    % The direction ofploting is clock-wise, so we change order of channels group
    reorderedCh = [FR CR TR PR OR OL PL TL CL FL];
    
    % we call the bottom nested function to draw each channel with its color on
    % the plot
    for i = 1:125
        draw1channel(i, Colors(infoCh(reorderedCh(i)),:), 125);
    end
    
    % Build points in 3D space for each channel
    points = zeros([length(MatrixData), 3]);
    for i = 1:length(MatrixData)
        points(i, 1) = 0;
        points(i, 2) = cos(t(i));
        points(i, 3) = sin(t(i));
    end
    
    all_n_conns = zeros([1, length(MatrixData)]);

%%  2. Plotting direct (positive) genaralization connectivity

    % First, check if there are direct (positive) generalization connections and
    % manage them
    if sum(sum(MatrixData > 0))

        % Prepare data to plot transforming in a 0 to 1 data weighted by
        % time length
        mask = MatrixData > 0; % !=
        MatrixMasked = abs(MatrixData) .* mask;
        Matrix = MatrixMasked(reorderedCh, reorderedCh)';
        Matrix = Matrix./lengthTime;
        
        
        % Now start lines
        % Find non-zero values of s and their indices
          [row,col,v] = find(Matrix);
          [v, indsV] =  sort(v, 'ascend');
          row = row(indsV);
          col = col(indsV);
          %%% Activate if you want white to color lines
          if gradient
          colormapwhite2black = [linspace(1,0.35, length(v))' linspace(1,0.35, length(v))' linspace(1,0.35, length(v))'];
          else 
          colormapwhite2black = [linspace(0.35,0.35, length(v))' linspace(0.35,0.35, length(v))' linspace(0.35,0.35, length(v))'];
          end  

        % Calculate line widths based on values of s (stored in v).
        minLineWidth  = 0.5; % minimum width of the line even if lowest decoding rate
        lineWidthCoef = 5; % coeficient of muliplier for width decoding rate
        lineWidth = v./max(v);
        if sum(lineWidth) == numel(lineWidth) && length(lineWidth) > 1 % all lines are the same width.
            lineWidth = repmat(minLineWidth,numel(lineWidth),1);
        else % lines of variable width depending on decoding rate.
            lineWidth = lineWidthCoef*lineWidth + minLineWidth;
        end
        
        
        % Create the trajectory of lines
        t      = (0.025: 0.05 :1)';
        t2  = [1-t, t].^2;
        
        clear Bx By
        Bx = zeros(length(v), 20);By = zeros(length(v), 20);
        % Loop over lines to be created to build them and plot them
        for idx = 1:length(v)
            % build trajectory
            Bx(idx, :)     = t2*[points(col(idx),2); points(row(idx),2)] ;
            By(idx, :)     = t2*[points(col(idx),3); points(row(idx),3)] ;
            % plot line
            plot3(ones(length(Bx(idx,:)),1)*0, Bx(idx,:),By(idx,:), 'LineWidth', lineWidth(idx), 'Color', colormapwhite2black(idx,:))
        end
    
        % Save number of connections for future calculation of size of each
        % channel dot
        n_conns = sum(Matrix > 0);
        all_n_conns = all_n_conns + n_conns;

    end  % End direct (positive) generalization connections

%%  3. Plotting inverse (negative) genaralization connectivity
    
    % Second, check if there are direct (positive) generalization connections and
    % manage them
    if sum(sum(MatrixData < 0))
        
        % Prepare data to plot transforming in a 0 to 1 data weighted by
        % time length
        mask = MatrixData < 0;
        MatrixMasked = abs(MatrixData) .* mask;
        Matrix = MatrixMasked(reorderedCh, reorderedCh)';
        Matrix = Matrix./lengthTime;
        
        % Now start lines
        % Find non-zero values of s and their indices
        [row,col,v] = find(Matrix);
        [v, indsV] =  sort(v, 'descend');
        row = row(indsV);
        col = col(indsV);
        %%% Activate if you want white to color lines 
        if gradient
        colormapwhite2red = [linspace(1,0.7, length(v))' linspace(1,0.35, length(v))' linspace(1,0.35, length(v))'];
        else
        colormapwhite2red = [linspace(0.7,0.7, length(v))' linspace(0.35,0.35, length(v))' linspace(0.35,0.35, length(v))'];
        end

        % Calculate line widths based on values of s (stored in v).
        minLineWidth  = 0.5; % minimum width of the line even if lowest decoding rate
        lineWidthCoef = 5; % coeficient of muliplier for width decoding rate
        lineWidth = v./max(v);
        if sum(lineWidth) == numel(lineWidth) && length(lineWidth) > 1 % all lines are the same width.
            lineWidth = repmat(minLineWidth,numel(lineWidth),1);
        else % lines of variable width depending on decoding rate.
            lineWidth = lineWidthCoef*lineWidth + minLineWidth;
        end
        
        
        % Create the trajectory of lines
        t      = (0.025: 0.05 :1)';
        t2  = [1-t, t].^2;
        
        clear Bx By
        Bx = zeros(length(v), 20);By = zeros(length(v), 20);
        % Loop over lines to be created to build them and plot them
        for idx = 1:length(v)
            % build trajectory
            Bx(idx, :)     = t2*[points(col(idx),2); points(row(idx),2)] ;
            By(idx, :)     = t2*[points(col(idx),3); points(row(idx),3)] ;
            % plot line
            plot3(ones(length(Bx(idx,:)),1)*0, Bx(idx,:),By(idx,:), 'LineWidth', lineWidth(idx), 'Color', colormapwhite2red(idx,:))
        end
    
        % Save number of connections for future calculation of size of each
        % channel dot
        n_conns = sum(Matrix > 0);
        all_n_conns = all_n_conns + n_conns;
        
    end  % End inverse (negative) generalization connections
    
    



%% 4. Draw each channel dot with the size for number of connections and color for decoding rate
    
    % Caluclate size depending on number of connections
    all_sizes = (all_n_conns/max(all_n_conns).*100) + 30;
    %scatter3(points(:,1), points(:,2), points(:,3), all_sizes, diag(MatrixData), 'filled'); clim([min(diag(MatrixData)), max(diag(MatrixData))]);
    % Draw points
    scatter3(points(:,1), points(:,2), points(:,3), all_sizes, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9]); 


%% 5. Visualize luminosity of figure, rotatation, labeling and legends

    
    light('Position',[-3 -2 -1]);
    light('Position',[3 1 -3]);
    axis vis3d; axis equal; view(3);
    set(gcf,'Color','w')
    view(-90,0)
    axis([-1.15 1.15 -1.15 1.15])
    axis off

    if legends

        % Brain Areas Labels
        text(0, 0, 1.2,'\bf Frontal','HorizontalAlignment','center','Color',Colors(1,:),'FontSize',12);
        c1 = text(0, -1.10, 0.55,'\bf R Central','HorizontalAlignment','center','Color',Colors(2,:),'FontSize',12); set(c1, 'Rotation', -60)
        t1 = text(0, -1.15, -0.4,'\bf R Temporal','HorizontalAlignment','center','Color',Colors(3,:),'FontSize',12); set(t1, 'Rotation', 70)
        p1 = text(0, -0.7, -1,'\bf R Parietal','HorizontalAlignment','center','Color',Colors(4,:),'FontSize',12); set(p1, 'Rotation', 30)
        text(0, 0, -1.2,'\bf Occipital','HorizontalAlignment','center','Color',Colors(5,:),'FontSize',12)
        c2 = text(0, 1.10, 0.55,'\bf L Central','HorizontalAlignment','center','Color',Colors(2,:),'FontSize',12); set(c2, 'Rotation', 60)
        t2 = text(0, 1.15, -0.4,'\bf L Temporal','HorizontalAlignment','center','Color',Colors(3,:),'FontSize',12); set(t2, 'Rotation', -70)
        p2 = text(0, 0.7, -1,'\bf L Parietal','HorizontalAlignment','center','Color',Colors(4,:),'FontSize',12); set(p2, 'Rotation', -30)
        
        % dots scatter
        scatter3(0, 1.15, -1.4, 30, [0.5712 0.8711 0.8889], 'filled')
        scatter3(0, 1.15, -1.55, 130, [0.5712 0.8711 0.8889], 'filled')


        % Legend Labels
        graytext = [0.2 0.2 0.2];
        text(0, 1.05, -1.4,"" + sprintf('%g', round(min(n_conns),2)) + " conns.",'HorizontalAlignment','left','Color',graytext,'FontSize',12)
        text(0, 1.05, -1.55,"" + sprintf('%g', round(max(n_conns),2)) + " conns.",'HorizontalAlignment','left','Color',graytext,'FontSize',12)
        
        %text(0, 0.10, -1.4,"" + sprintf('%g', round(t,2)) + " time",'HorizontalAlignment','left','Color',graytext,'FontSize',12)
        text(0, 0.10, -1.4,"Time:",'HorizontalAlignment','left','Color',graytext,'FontSize',12)
        text(0, 0, -1.55,"\bf " + sprintf('%g', round(timestep)) + " to " + sprintf('%g', round(timestep + lengthTime - 1)),'HorizontalAlignment','center','Color',graytext,'FontSize',14);
    
    end


%% 6. Auxiliar function for creating the visualization of channels groups in outer ring
    
    
    function draw1channel(numCh, colorCh, totalCh)
        
        t = linspace(pi/2,2*pi + pi/2, totalCh);
        rin = 1.05; % radio interior
        rout = 1.12; % radio exterior
        center = [0, 0]; 
        xin = rin*cos(t);
        xout = rout*cos(t);
        yin = rin*sin(t);
        yout = rout*sin(t);
        z1 = 0-0.035;
        z2 = 0+0.035;

        hold on;
    
        % check if is not the last part
        if numCh ~= totalCh  % Only changes that each channel goes from 1 to 2 position and in last one, from end to first
            % Draw bottom and top parts of the block
            bottom = patch(z1*ones(1,4), ...
                           center(1)+[xin(numCh:numCh+1) xout(numCh+1) xout(numCh)], ...
                           center(2)+[yin(numCh:numCh+1) yout(numCh+1) yout(numCh)], ...
                           colorCh);
            top = patch(z2*ones(1,4), ...
                           center(1)+[xin(numCh:numCh+1) xout(numCh+1) xout(numCh)], ...
                           center(2)+[yin(numCh:numCh+1) yout(numCh+1) yout(numCh)], ...
                           colorCh);
       
            % Create Cylinder
            X = [cos(t); cos(t)];
            Y = [sin(t); sin(t)];
            Z = [zeros([1,length(t)+1]); ones([1,length(t)+1])];
        
        
            % Draw outer part of block
            tmp2 = rout*X+center(1);
            tmp3 = rout*Y+center(2);
            tmp1 = Z*(z2-z1)+z1;
            outer = surf(tmp1(:,numCh:numCh+1), ...
                         tmp2(:,numCh:numCh+1), ...
                         tmp3(:,numCh:numCh+1));
            
            % Draw inner part of block
            tmp2 = rin*X+center(1);
            tmp3 = rin*Y+center(2);
            tmp1 = Z*(z2-z1)+z1;
            inner = surf(tmp1(:,numCh:numCh+1), ...
                         tmp2(:,numCh:numCh+1), ...
                         tmp3(:,numCh:numCh+1));
            
        else  % Only changes that each channel goes from 1 to 2 position and in last one, from end to first
            % Draw bottom and top parts of the block
            bottom = patch(z1*ones(1,4), ...
                       center(1)+[xin(numCh) xin(1) xout(1) xout(numCh)], ...
                       center(2)+[yin(numCh) yin(1) yout(1) yout(numCh)], ...
                       colorCh);
            top = patch(z2*ones(1,4), ...
                       center(1)+[xin(numCh) xin(1) xout(1) xout(numCh)], ...
                       center(2)+[yin(numCh) yin(1) yout(1) yout(numCh)], ...
                       colorCh);

            % Draw outer part of block
            [X,Y,Z] = cylinder(1,length(xin));
            tmp2 = rout*X+center(1);
            tmp3 = rout*Y+center(2);
            tmp1 = Z*(z2-z1)+z1;
            outer = surf([tmp1(:,numCh) tmp1(:,1)], ...%rout*X+center(1), ...
                         [tmp2(:,numCh) tmp2(:,1)], ...%rout*Y+center(2), ...
                         [tmp3(:,numCh) tmp3(:,1)]);%Z*(z2-z1)+z1);

            % Draw inner part of block
            tmp2 = rin*X+center(1);
            tmp3 = rin*Y+center(2);
            tmp1 = Z*(z2-z1)+z1;
            inner = surf(tmp1(:,numCh:1), ...
                         tmp2(:,numCh:1), ...
                         tmp3(:,numCh:1));
    
        end
        
        % Set color for the channel block
        set([bottom, top, outer, inner], ... 
            'FaceColor', colorCh, ...
            'FaceAlpha', 0.99, ...
            'linestyle', 'none', ...
            'SpecularStrength', 0.7);
    end
    

end
