close all;
clf;

global array
global arrsize
global chemotherapy
global stop
global TIME
global chemotherapy_dead

arrsize=100; %100x100 array 
chemotherapy=0;
stop=0;
TIME=1;
chemotherapy_dead=0;

%creates the array
for ii=1:arrsize
    for jj=1:arrsize 
        array(ii,jj,1)=0; %cellular state
        %0=empty, 1=proliferative, 2=necrotic, 3=artery, 4=hypoxic,
        %5=capillary, 6=Chemo Necrotic
        array(ii,jj,2)=0; %VEGF
        array(ii,jj,3)=0.75; %oxygen
        array(ii,jj,4)=0.75; %nutrients
        array(ii,jj,5)=0; %chemotherapy drug
        array(ii,jj,6)=0; %age of cell
        array(ii,jj,7)=0; %drug resistance state
        array(ii,jj,8)=0; %VEGF old
        array(ii,jj,9)=0.75; %oxygen old
        array(ii,jj,10)=0.75; %nutrients old
        array(ii,jj,11)=0; %chemotherapy drug old
        array(ii,jj,12)=0; %cellular state old
        array(ii,jj,13)=0; %hypoxia "age"
        array(ii,jj,14)=0; %KRAS mutation state
    end
end

ind=floor(arrsize*.1/2);
ind=arrsize/2-ind:arrsize/2+ind; %finding middle of lattice
for ii=ind
    for jj=ind
        if(rand(1)>.1) % lattice points within the square have a 90% probability of being filled with a cancer cell
            array(ii,jj,1)=1;
            array(ii,jj,12)=1;
            array(ii,jj,3)=.5; %occupied cells have decreased oxygen and nutrients to start
            array(ii,jj,4)=.5;
            array(ii,jj,9)=.5;
            array(ii,jj,10)=.5;
        end
    end
end

array((10:16:90),2:3,1)=5; %making capillaries
array((10:16:90),2:3,12)=5;
array((10:16:90),arrsize-2:arrsize-1,1)=5;
array((10:16:90),arrsize-2:arrsize-1,12)=5; 
array((10:80:90),4:6,1)=5;
array((10:80:90),4:6,12)=5;
array((10:80:90),arrsize-5:arrsize-3,1)=5;
array((10:80:90),arrsize-5:arrsize-3,12)=5; 
array(2:3,(20:20:80),1)=5;
array(2:3,(20:20:80),12)=5;
array(arrsize-2:arrsize-1,(20:20:80),1)=5;
array(arrsize-2:arrsize-1,(20:20:80),12)=5; 
array(1,:,1)=3;%making arteries
array(1,:,12)=3;
array(arrsize,:,1)=3;
array(arrsize,:,12)=3;
array(:,1,1)=3;
array(:,1,12)=3;
array(:,arrsize,1)=3;
array(:,arrsize,12)=3;


proliferative=zeros(1,100);
hypoxic=zeros(1,100);
necrotic=zeros(1,100);
capillary_endothelial=zeros(1,100);
chemotherapy_resistant=zeros(1,100);
tumor_diameter=zeros(1,100);
necrotic_diameter=zeros(1,100);
chemo_dead=zeros(1,100);

for t = 1:100
    live();
    [row,col] = find(array(:,:,1) == 1); %finding row and column index for every alive cell
    idx = find(array(:,:,1) == 1); % finding matrix index of alive cell locations
    [ii,jj] = ind2sub(size(array),idx); %getting x and y coordinates of proliferative cells
    [row,col] = find(array(:,:,1) == 2); 
    idx = find(array(:,:,1) == 2);
    [kk,ll] = ind2sub(size(array),idx);
    [row,col] = find(array(:,:,1) == 3);
    idx = find(array(:,:,1) == 3);
    [mm,nn] = ind2sub(size(array),idx);
    [row,col] = find(array(:,:,1) == 4);
    idx = find(array(:,:,1) == 4);
    [bb,cc] = ind2sub(size(array),idx);
    [row,col] = find(array(:,:,1) == 5);
    idx = find(array(:,:,1) == 5);
    [dd,ee] = ind2sub(size(array),idx);
    idx = find(array(:,:,1) == 6);
    [qq,uu] = ind2sub(size(array),idx);
    figure(1); clf();
    hold on
    scatter(ii,jj,'o','MarkerFaceColor','g','MarkerEdgeColor','g'); %green cells are proliferative
    scatter(kk,ll,'o','MarkerFaceColor','k','MarkerEdgeColor','k'); %necrotic cells are black
    scatter(mm,nn,'o','MarkerFaceColor','k','MarkerEdgeColor','r'); %artery endothelial cells are black with red outline
    scatter(bb,cc,'o','MarkerFaceColor','b','MarkerEdgeColor','b'); %hypoxic cells are blue
    scatter(dd,ee,'o','MarkerFaceColor','r','MarkerEdgeColor','b'); %capillary endothelial cells are red with red outline
    scatter(qq,uu,'o','MarkerFaceColor',[1 0.75 0],'MarkerEdgeColor','b'); %capillary endothelial cells are red with red outline
    axis([0 100 0 100])
    hold off
    drawnow
    proliferative(t)=sum(sum(array(:,:,1)==1));
    hypoxic(t)=sum(sum(array(:,:,1)==4));
    necrotic(t)=sum(sum(array(:,:,1)==2));
    capillary_endothelial(t)=sum(sum(array(:,:,1)==5));
    chemotherapy_resistant(t)=sum(sum(array(:,:,7)>=1));
    chemo_dead(t)=chemotherapy_dead;
    a=find(array(:,arrsize/2,1)==1 | array(:,arrsize/2,1)==4 | array(:, arrsize/2,1)==2); %finding proliferative tumor cells that occupy the vertical midline of the grid; this will be used to quantify tumor diameter
    if (max(a)-min(a)) > 0
        tumor_diameter(t)=max(a)-min(a);
    else
        tumor_diameter(t)=0;
    end
    b=find(array(:,arrsize/2,1)==2); %finding necrotic tumor cells that occupy the vertical midline of the grid; this will be used to quantify necrotic core diameter
    if (max(b)-min(b)) > 0
        necrotic_diameter(t)=max(b)-min(b);
    else
        necrotic_diameter(t)=0;
    end
    if rem(t,5)==0 %every fifth time point, plot gradients of VEGF, oxygen, and nutrients; darker color=higher concentration
        VEGF_plot=array(:,:,2)/(max(max(array(:,:,2))));
        figure;
        for ff=0:9
            [row,col] = find(VEGF_plot >=ff/10 & VEGF_plot <= (ff/10+.1)); 
            idx = find(VEGF_plot >=ff/10 & VEGF_plot <= (ff/10+.1)); 
            [gg,ee] = ind2sub(size(VEGF_plot),idx); 
            title(['VEGF Gradient after ',num2str(t),' Time Points']);
            hold on
            scatter(gg,ee,'s','MarkerFaceColor',[0 0 abs(ff-9)/10],'MarkerEdgeColor',[0 0 abs(ff-9)/10]);
            axis([0 100 0 100])
            drawnow
        end
        hold off
        oxygen_plot=array(:,:,3)/(max(max(array(:,:,3))));
        figure;
        for ff=0:9
            [row,col] = find(oxygen_plot >=ff/10 & oxygen_plot <= (ff/10+.1)); 
            idx = find(oxygen_plot >=ff/10 & oxygen_plot <= (ff/10+.1)); 
            [gg,ee] = ind2sub(size(oxygen_plot),idx); 
            title(['Oxygen Gradient after ',num2str(t),' Time Points']);
            hold on
            scatter(gg,ee,'s','MarkerFaceColor',[abs(ff-9)/10 0 0],'MarkerEdgeColor',[abs(ff-9)/10 0 0]);
            axis([0 100 0 100])
            drawnow
        end
        hold off
        chemotherapy_drug_plot=array(:,:,5)/(max(max(array(:,:,5))));
        figure;
        for ff=0:9
            [row,col] = find(chemotherapy_drug_plot >=ff/10 & chemotherapy_drug_plot <= (ff/+.1)); 
            idx = find(chemotherapy_drug_plot >=ff/10 & chemotherapy_drug_plot <= (ff/10+.1)); 
            [gg,ee] = ind2sub(size(chemotherapy_drug_plot),idx); 
            title(['Chemotherapy Drug Gradient after ',num2str(t),' Time Points']);
            hold on
            scatter(gg,ee,'s','MarkerFaceColor',[0 abs(ff-9)/10 0],'MarkerEdgeColor',[0 abs(ff-9)/10 0]);
            axis([0 100 0 100])
            drawnow
        end
        hold off
    end    
    if(sum(array(:,1,1)==1)>0||sum(array(:,arrsize,1)==1)>0||sum(array(1,:,1)==1)>0||sum(array(arrsize,:,1)==1)>0)
        break %if a proliferative tumor cell touches an outer boundary, end the simulation
    end
    if stop==1 %if an artery cell touches a proliferative tumor cell, end the simulation
        break
    end
    if sum(sum(array(:,:,1)==1)) ==0 && sum(sum(array(:,:,1)==4))==0 %end the simulation if all cells are necrotic
        break
    end
    if chemotherapy==0 %when tumor diameter is greater than 50 lattice points (2.5 cm), start chemotherapy treatment
       if sum(sum(array(:,:,1)==1))+sum(sum(array(:,:,1)==2))+sum(sum(array(:,:,1)==4)) > 2000
           beg=TIME;
           chemotherapy=1
       end
    else
        if(beg-TIME>30|sum(sum(array(:,:,1)==1))+sum(sum(array(:,:,1)==4)) < 50)
            chemotherapy
        end
    end
    t
    if t<100
        TIME=TIME+1;
    end
end

total_cells=hypoxic+proliferative+necrotic;
necrotic_fraction=(necrotic(1:TIME))./(total_cells(1:TIME));
viable_fraction=(proliferative(1:TIME)+hypoxic(1:TIME))./(total_cells(1:TIME));
necrotic_diameter=necrotic_diameter.*25;
tumor_diameter=tumor_diameter.*25;

figure;
plot(1:TIME,total_cells(1:TIME),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Tumor Cell Number vs. Time');
xlabel('Time Steps');
ylabel('Tumor Cell Number');

figure;
plot(1:TIME,capillary_endothelial(1:TIME),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Capillary Endothelial Cells vs. Time');
xlabel('Time Steps');
ylabel('Capillary Endothelial Cell Number');

figure;
plot(1:TIME,necrotic_fraction, 'o','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on;
plot(1:TIME,viable_fraction,'o','MarkerFaceColor','g','MarkerEdgeColor','g');
title('Necrotic and Viable Fractions vs. Time');
xlabel('Time Steps');
ylabel('Percentage of Total Tumor Cells');
legend('Necrotic', 'Viable');
hold off;

figure;
plot(1:TIME,tumor_diameter(1:TIME),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Tumor Diameter vs. Time');
xlabel('Time Steps');
ylabel('Tumor Diameter (?m)');


figure;
plot(1:TIME,necrotic_diameter(1:TIME),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Necrotic Diameter vs. Time');
xlabel('Time Steps');
ylabel('Necrotic Diameter (?m)');


figure;
plot(1:TIME,chemotherapy_resistant(1:TIME),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Chemotherapy Resistant Cells vs. Time');
xlabel('Time Steps');
ylabel('Number of Chemotherapy-Resistant Cells');


figure;
plot(1:TIME,chemo_dead(1:TIME),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Cells Killed by Chemotherapy vs. Time');
xlabel('Time Steps');
ylabel('Number of Cells Killed by Chemotherapy');


%Function to run on each cell to do all functions
function live()
    global arrsize
    global array
    for ii=1:arrsize
        for jj=1:arrsize
            diffusion(ii,jj);
            if(array(ii,jj,1)==0) %if the lattice point is empty, run grow
                grow(ii,jj);
            end
            consume(ii,jj);
        end
    end
    array(:,:,8)=array(:,:,2); %updating VEGF for next time step
    array(:,:,9)=array(:,:,3); %oxygen
    array(:,:,10)=array(:,:,4); %nutrients 
    array(:,:,11)=array(:,:,5); %drugs
    array(:,:,12)=array(:,:,1); %cellular state
end

function grow(x,y) %populate empty lattices with new cell
    global array
    [neighcells, ~]=neighbors(x,y); %calculating neighboring cancer cells and endothelial cells;neighcells/vessels is 6.83 at max
    if (neighcells>=1) && (array(x,y,3)>0.18) && (array(x,y,4)>.18) && (rand(1)>0.2) && array(x,y,5)<0.25 %populates an empty lattice with a cell under these conditions: at least one vertical/horizontal neighbor (or two diagonal neighbors), sufficient oxygen and nutrients
        array(x,y,1)=1;
    end
end

function growves(x,y,time) %angiogenesis function
    global arrsize
    global array
    persistent log
    count=1;
    [~, neighvessels]=neighbors(x,y); %calculating number of neighboring vessels
    if array(x,y,12)==5 && time ~= 0 %if the lattice point is a capillary endothelial cell
        if(rand(1)>.999) && array(x,y,6)>10 %branching can occur if capillary endothelial cells are more than 10 time steps old
            n=1+(1/sqrt(1));
        else
            n=1;
        end
        if array(x,y,2)>0.25 && rand(1)>0.1 && neighvessels<=n && neighvessels>=(1/sqrt(2)) %If VEGF is sufficiently large and the endothelial cell is at the tip of a blood vessel (i.e. has one diagonal or one horizontal neighbor)
            for ii=(x-1):(x+1)
                for jj=(y-1):(y+1)
                    if rangeBox(ii,jj) %survery VEGF concentrations in the Moore neighborhood around the endothelial cell
                        [~, neighvessels]=neighbors(ii,jj);
                        if array(ii,jj,12)==3 || array(ii,jj,12)==5 || array(ii,jj,1)==1 || array(ii,jj,1)==2 || array(ii,jj,1)==4 || neighvessels > 1 %don't look at VEGF concentration if the cell is already occupied or if it has lots of capillary endothelial cells neighbors
                            VEGF(count)=0;
                            count=count+1;
                        else
                            VEGF(count)=array(ii,jj,2); %store VEGF concentration for this lattice point
                            count=count+1;
                        end
                    else
                        VEGF(count)=0;
                        count=count+1;
                    end
                end
            end
            maximumVEGF=max(VEGF); %find max VEGF concentration in the Moore neighborhood
            counter=0;
            if maximumVEGF > 0
                for ii=(x-1):(x+1)
                    for jj=(y-1):(y+1)
                        if rangeBox(ii,jj) && counter==0
                            if array(ii,jj,1)==0 %if the cell is empty
                                if array(ii,jj,2)==maximumVEGF %and is the cell with the max VEGF concentration in the neighborhood
                                    if(isempty(log))
                                        log=[];
                                    end
                                    if(sum(log==sub2ind(size(array),x,y))==0 || time~=3)
                                        array(ii,jj,1)=5; %the lattice point with the maximum VEGF concentration becomes an endothelial cell
                                        array(ii,jj,12)=5;
                                        log=[log sub2ind(size(array),ii,jj)];
                                        growves(ii,jj,time-1);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if(x==arrsize&&y==arrsize)
        log=[];
    end
end

function consume(x,y) %consumption and production of diffusable factors, hypoxia and death induction
    global arrsize
    global array
    global chemotherapy
    global stop
    global chemotherapy_dead
    if(array(x,y,1)==6 & array(x,y,6) > 2)
        array(x,y,1)=0;
    end
    growves(x,y,3); %perform angiogenesis
    if(array(x,y,12)==3) %if the lattice point is an artery
       [neighcells, ~]=neighbors(x,y); %if the artery has any proliferative cells as its neighbor
       if neighcells > 0
           stop=1;
       end
       array(x,y,2)=0; %VEGF at that lattice point is consumed by artery
       if rand(1)>0.5 %only produce factors 50% of the time
           for ii=x-10:x+10
               for jj=x-10:x+10
                   if(rangeBox(ii,jj))
                       array(x,y,3)=2; %Produces oxygen
                       array(x,y,4)=2; %Produces nutrients
                       if chemotherapy==1
                           array(x,y,5)=4; %Produces drug when chemotherapy is activated
                       end
                   end
               end
           end
        end
    end
    if(array(x,y,12)==5) %if the lattice point is a capillary, it has greater perfusion compared to arteries
        array(x,y,2)=0; %VEGF at that lattice point is consumed by artery
        if (array(x,y,6)>2) %if the capillary is not "immature", which we define as being 2 or less timesteps old
            if rand(1)>0.5 %only produce factors 50% of the time
                for ii=x-10:x+10
                    for jj=x-10:x+10
                        if(rangeBox(ii,jj))
                            array(x,y,3)=5; %Produces oxygen
                            array(x,y,4)=5; %Produces nutrients
                            if chemotherapy==1
                                array(x,y,5)=7;
                            end
                        end
                    end
                end
            end
        else %if the capillary cell is immature
            if rand(1)>0.5 %only produce factors 50% of the time
                for ii=x-10:x+10
                    for jj=x-10:x+10
                        if(rangeBox(ii,jj))
                            array(x,y,3)=2; %Produces oxygen
                            array(x,y,4)=2; %Produces nutrients
                            if chemotherapy==1
                                array(x,y,5)=4;
                            end
                        end
                    end
                end
            end            
        end
     end
    array(x,y,2)=array(x,y,2)*.95; %Degradtion of VEGF
    array(x,y,5)=array(x,y,5)*.9; %Degradation of drug
    array(x,y,3)=array(x,y,3)*.975; %Degradation of oxygen: higher b/c other, healthy cells in the area are consuming it
    array(x,y,4)=array(x,y,4)*.975; %Same logic for nutrients
    if(array(x,y,1)==1 || array(x,y,1)==4) %if the point is occupied by a proliferative or hypoxic cell
       ox=array(x,y,3)-neighborsnew(x,y)/25; %degrade oxygen and nutrients based on number of neighbors
       if ox < 0 %don't let this fall below 0
           ox=0;
       end
       nut=array(x,y,4)-neighborsnew(x,y)/25;
       if nut < 0
           nut=0;
       end
       if(ox<.18) %if the cell has low oxygen levels
           array(x,y,1)=4; %label it hypoxic
           if array(x,y,12)==4
               array(x,y,13)=array(x,y,13)+1.25+(.5-ox)+(nut/5); %if it was also hypoxic at the last time point, increase its hypoxic age according to oxygen and nutrient levels, and by a constant factor (1.25) each time step
           else
               array(x,y,13)=1; %if it's newly hypoxic, assign its hypoxic age as 1
           end
           if rand(1)>(1-(array(x,y,13)/100)-(array(x,y,6)/100)) 
               array(x,y,14)=1; %the hypoxic cell is hit with a KRAS mutation allowing it to pump out more VEGF, and a rate consistent with the amount of time its been hypoxic and its overall cellular age
           end
           if rand(1)>0.1 && array(x,y,14)==1 %if the cell is KRAS mutated
               for ii=x-15:x+15  %cycling through the lattice points surrounding the middle point
                   for jj=y-15:y+15 
                       if(rangeBox(ii,jj))
                           array(ii,jj,2)=array(ii,jj,2)+15; %Produce more VEGF than normal
                       end
                   end
               end
           else %if the cell is not KRAS mutated
               for ii=x-15:x+15  %cycling through the lattice points surrounding the middle point
                   for jj=y-15:y+15 
                       if(rangeBox(ii,jj))
                           array(ii,jj,2)=array(ii,jj,2)+10; %Produce normal VEGF 
                       end
                   end
               end
           end              
       else
           array(x,y,1)=1; %if oxygen is not low, cell is proliferative and its hypoxic timer is reset to 0
           array(x,y,13)=0;
       end
    end
    if array(x,y,13) >= 23 %if hypoxic age of the cell is too high, the cell becomes necrotic 
        array(x,y,1)=2;
    end
    if (array(x,y,1)==1 || array(x,y,1)==4) && array(x,y,5) > 0.25 %cell dies if chemotherapy concentration in its lattice point is too high
        array(x,y,5)=array(x,y,5)-neighborsnew(x,y); %consumption of drug by the cell
        if rand(1)>0.12 && array(x,y,7)<1 | rand(1)>0.3 && array(x,y,1)>=1 %cell has 88% probability of dying and being removed from the simulation
            array(x,y,1)=6;
            array(x,y,6)=0;
            array(x,y,7)=0;
            array(x,y,13)=0;
            array(x,y,14)=0;
            chemotherapy_dead=chemotherapy_dead+1;
        else %if cell survives the treatment
            if array(x,y,7) < 1
                array(x,y,7)=array(x,y,7)+(array(x,y,6)/70); %its will develop chemotherapeautic resistance if it survives multiple treatments; less treatments need to be survived if the cell is older in age, reflecting the accumulation of mutations over time
                if array(x,y,7) > 1
                    array(x,y,7)=1;
                end
            end
        end
    end
    if(array(x,y,1)~=0) 
          array(x,y,6)=array(x,y,6)+1; %if the non-necrotic tumor cell or blood vessel cell doesn't die, increment its age
    end
    array(x,y,3)=array(x,y,3)-neighborsnew(x,y)/25; %Reduce oxygen depending on the number of alive neighbors
    if array(x,y,3) < 0
        array(x,y,3)=0;
    end
    array(x,y,4)=array(x,y,4)-neighborsnew(x,y)/25; %Reduce nutrients depending on the number of alive neighbors
    if array(x,y,4) < 0
        array(x,y,4)=0;
    end
end
       
       %Calculates the number of neighbors
function [cancer, vessel] = neighbors(x, y)
    global array
    vessel=0;
    cancer=0;
    for ii=(x-1):(x+1)
        for jj=(y-1):(y+1) %Moore neighborhood: the 8 cells surrounding the lattice point of interest
            if(rangeBox(ii,jj)) %makes sure you don't calculate neighbors "outside" of the lattice grid
                if~(ii==x && jj==y) %excludes the middle cell from the calculation
                   if(array(ii,jj,12)==1) %if a proliferative cell exists in the surrounding position; hypoxic cells are considered to be quiscent and don't contribute to the neighbor count
                       cancer=cancer+(1/(distanceeqn(x,y,ii,jj)*(10^6))); %vertical/horizontal neighbors contribute 1; diagonal neighbors contribute 1/sqrt(2)
                   elseif(array(ii,jj,12)==5) %if a capillary endothelial cell exists in the surrounding position
                       vessel=vessel+(1/(distanceeqn(x,y,ii,jj)*(10^6)));
                   end
                end
            end
        end
    end
end

function cancer = neighborsnew(x, y)
    global array
    cancer=0;
    for ii=x-1:x+1 
        for jj=y-1:y+1 
            if(rangeBox(ii,jj)) %makes sure you don't calculate neighbors "outside" of the lattice grid
                if~(ii==x && jj==y) %excludes the middle cell from the calculation
                   if(array(ii,jj,1)==1) %if a proliferative cell exists in the surrounding position
                       cancer=cancer+(1/(distanceeqn(x,y,ii,jj)*(10^6)));
                   end
                end
            end
        end
    end
end

% Csalculation of the diffusion for point (x,y)
% Only the (x,y) lattice point is updated; all others around it are not,
% until their lattice point is being analyzed specifically
function diffusion(x,y)
    global array
    for ii=x-3:x+3 %cycling through lattice points surrounding the middle point
        for jj=y-3:y+3 
            if(rangeBox(ii,jj))
                if~(ii==x && jj==y) %excluding the middle point itsel
                    delt=diffuse(x,y,ii,jj,8); %calculated diffusion using Fick's Law
                    array(x,y,2)=array(x,y,2)+delt; %Value is updated in the new matrix
                    if array(x,y,2) < 0 %don't let values go below 0
                        array(x,y,2)=0;
                    end
                    delt=diffuse(x,y,ii,jj,9);
                    array(x,y,3)=array(x,y,3)+delt;
                    if array(x,y,3) < 0
                        array(x,y,3)=0;
                    end
                    delt=diffuse(x,y,ii,jj,10);
                    array(x,y,4)=array(x,y,4)+delt;
                    if array(x,y,4) < 0
                        array(x,y,4)=0;
                    end
                    delt=diffuse(x,y,ii,jj,11);
                    array(x,y,5)=array(x,y,5)+delt;
                    if array(x,y,5) < 0
                        array(x,y,5)=0;
                    end
                end
            end
        end
    end
end

%makes sure the lattice point you're evaluating is in the array
function inside = rangeBox(x,y)
    global arrsize
    inside=false;
    if(x>=1 && x<=arrsize)
            if(y>=1 && y<=arrsize)
                inside=true;
            end
    end
end

%calculates the eculidian distance between two points
function dist = distanceeqn(x,y,x2,y2)%x and y is the lattice point you're at; x2 and y2 is the lattice point in that lattice's neighborhood that you're evaluating
    dist=sqrt((x-x2)^2+(y-y2)^2)*(10^-6); %distance is either 1 or sqrt(2); in micrometers
end

%calculates the flux of a diffusable factor between 2 points
function flux = diffuse(x,y,x2,y2,type) %x and y is the lattice point you're at; x2 and y2 is the lattice point in that lattice's neighborhood that you're calculating diffusion for
    global array
    if(type==8||type==9||type==10||type==11)
        if (type==8)
            permeability=1*(10^-1); %Permeability of VEGF (^-11)
        elseif (type==9)
            permeability=2.5*(10^-2); %Permeability of oxygen in tissue; m^2/s (^-9)
        elseif (type==10)
            permeability=1.63*(10^-2); %Permeability of glucose (^-10)
        else
            permeability=1*(10^-2); %Permeability of hypoxia-activated chemotherapy drugs (^-10)
        end
        flux=-permeability*((array(x,y,type)-array(x2,y2,type))/(distanceeqn(x,y,x2,y2)*10^6)); %Fick's Law of Diffusion
    else
        flux=0;
    end
end