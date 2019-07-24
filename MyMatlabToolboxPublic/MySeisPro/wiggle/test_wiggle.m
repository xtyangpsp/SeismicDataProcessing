%% Basic data
y = 200 + 4*(1:200);
x = 100 + 25*(1:20);
C = 2*randn(length(y),length(x));
D = filter([1,1,1],1,C);

%% Testing axis reescaling
figure(1); clf
subplot(131); wiggle(D);              title('Default');                  
subplot(132); wiggle(y, x, D);        title('Reescaled axes');                
subplot(133); wiggle(y, x, D, '-');   title('Negative lobes filled'); 

%% Testing different colors for lines and lobes
figure(2); clf
subplot(131); wiggle(y, x, D);        title('Default');                  
subplot(132); wiggle(y, x, D, 'r');   title('Red lines');                
subplot(133); wiggle(y, x, D, 'rM');  title('Red lines, magenta lobes'); 

%% Testing wiggle direction (vertical x horizontal)
figure(3); clf
subplot(131); wiggle(y, x, D);        title('Default');     
subplot(132); wiggle(y, x, D, '2');   title('2x extrusion');    
subplot(133); wiggle(y, x, D, '2h');  title('2x extrusion and horiz. lines'); 

%% Testing more complex combinations
figure(4); clf
subplot(131); wiggle(y, x, D, 'iBR'); title('Blue and red lobes and no lines'); 
subplot(132); wiggle(y, x, D, '-Bc'); title('Blue lobes on the left, cyan lines');     
subplot(133); wiggle(y, x, D, 'rI');  title('Red lines, no lobes'); 

%% Testing wiggles over image
figure(5); clf
subplot(131); wiggle(y, x, D);                                title('Default');                  
subplot(132); imagesc(x,y,D); hold on; wiggle(y, x, D);       title('Default wiggles over image');             
subplot(133); imagesc(x,y,D); hold on; wiggle(y, x, D, 'I');  title('Black lines only over image');  
colormap bwr

%% Testing colormapped lobes and extrusion
figure(6); clf
subplot(131); wiggle(y, x, D, 'BR');      title('Blue and red lobes');                  
subplot(132); wiggle(y, x, D, '1.25**');  title('Black lines with colormapped lobes'); 
subplot(133); wiggle(y, x, D, '1.5XX');  title('Black lines only over image');                
colormap bwr

