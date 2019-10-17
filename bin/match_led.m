function [real_sphere,real_status,rearranged] = match_led(vir_sphere,vir_status,lid)
    % Author: Yue Zhang   Date: 02/01/2018
    % Note: 
    %   ALL SUB-FUNCTIONS in this script (everything except "match_led")
    %   are either MODIFIED OR COPIED from the script "sphereEachLed.m"
    %   written by Florian A.Dehmelt.
    %
    % The major usage of this function is to find the corresponding pixel
    % values for each LED on the spherical arena on the virtual sphere
    % Input:
    %      vir_sphere: 
    %           M-by-3 matrix, the location of the M points on 
    %           the virtual sphere.
    %      vir_status: 
    %           M-by-N matrix, the on/off status of each point 
    %           on the virtual sphere at all N frames
    %      lid: 
    %           can be 'none','top only','bottom only', used to define the
    %           structure of the spherical arena (presence and position
    %           of the lid)
    % Output:
    %       real_sphere:
    %           M-by-3 matrix, the location of the M points on the real
    %           sphere, if input var 'lid' == 'none', M = 15104;
    %       real_status:
    %           M-by-N matrix, the on/off status for each LED on the arena
    %       rearranged:
    %           64x240xN matrix, the matrix can be saved in CF card for
    %           programing the arena
    %

    if ~exist('lid','var'), lid = 'none'; end
    
    real_sphere = sphereLed(lid); 
    real_shape  = size(real_sphere);
    real_shape  = [real_shape(2:end),size(vir_status,2)];
    real_sphere = reshape(real_sphere,3,[]);
    real_radius = mean(sqrt(sum(real_sphere.^2,1)));
    norm_sphere = (real_sphere )./real_radius;
    vir_order = knnsearch(vir_sphere,norm_sphere');
    real_status = vir_status(vir_order,:);
    real_status_r = permute(reshape(real_status,real_shape),[1,2,4,3]);
    rearranged = sphereFakeCylinder(real_status_r,lid);
    real_sphere = real_sphere';
end


function ledposcartesian = sphereLed(lid)
% Modified from Florian's "sphereEachLed" function
% Used to compute the location of each LED on the arena
if ~exist('lid','var'), lid = 'none'; end
switch lid
    case 'none'
      numtile = [5 9 11 13 14 15 14 13 11 8 5];
    case 'top only'
      numtile = [2 5 9 11 13 14 15 14 13 11 8 5];
    case 'bottom only'
      numtile = [5 9 11 13 14 15 14 13 11 8 5 2];
    otherwise
      error(['-- How many lids are covering the poles of the sphere? ', ...
             'Three scenarios are possible: ''none'', ''top only'', ', ...
             'and ''bottom only''. Please select one of them. --']);
end

[sphereradius, ribbonradius, ribbonangle, ribangle] = sphereShape(lid);

tilepos = sphereTilePosition(ribbonradius, ribbonangle, ribangle, ...
                           numtile, lid);

[ledposcartesian,~] = sphereLedPosition(tilepos,sphereradius);


ledposcartesian2 = ledposcartesian .* repmat([1;-1;1], [1, ...
                                   size(ledposcartesian,2), ...
                                   size(ledposcartesian,3), ...
                                   size(ledposcartesian,4)]);
ledposcartesian2 = flipdim(ledposcartesian2,2);

ledposcartesian = cat(4, ledposcartesian, ledposcartesian2);                       
end

function rearranged = sphereFakeCylinder(original,lid)
  
  % FIRST, REORDER TILES BY TILE ID NUMBER
  %
  % Rearrange information on the activity of each individual LED based on
  % the ID number of the tile upon which they are located. These tile
  % numbers (ranging from 1 to 240) are displayed on the spherical arena.

  % List the ID numbers of all tiles in one hemisphere, top-to-bottom and
  % meridian-to-side (i.e., going through one horiz. row after another).
  hemisphere1 = ...
    [115, 116,  65,  45,  47, ...
     113, 114,  80,  78,  66,  46, 48,  4,  2, ...
     119, 120, 118,  79,  77,  67, 56, 15, 14, 16,  3, ...
     112, 111, 110, 109, 117,  61, 68, 54, 38, 37, 12, 13,  1, ...
     108, 107, 106, 105,  85,  86, 63, 52, 55, 42, 39, 11, 24,  6, ...
      84,  83,  82,  81,  88,  87, 64, 51, 53, 43, 40,  9, 10, 23, 5, ...
     104, 103, 102, 101,  74,  76, 62, 50, 44, 41, 25, 18, 17, 22, ...
     100,  98,  97,  90,  73,  75, 49, 28, 27, 26, 20, 19, 21, ...
      99,  92,  89,  72,  71,  60, 58, 33, 34, 30, ...
      29,  94,  91,  70,  69,  59, 36, 35, 31, ...
      95,  96,  93,  57,  32];
  
  % List the ID numbers of the tiles on the lid of the same hemisphere.
  lid1 = [7, 8];
  
  
%   % The following are dummy IDs created for the second hemisphere and the
%   % second half of the lid. They must be replaced with the true IDs there
%   % as soon as Kun Wang has made these available.
%   hemisphere2 = hemisphere1 + 120;
%   lid2        = lid1        + 120;

  hemisphere2 = 120 + ...
    [100, 71, 72, 32, 28, ...
      99, 69, 70, 29, 30, 31, 26,  6,  8, ...
      98, 97, 88, 87, 68, 22, 23, 24, 25, 27,  7, ...
      91, 92, 85, 86, 66, 67, 21, 19, 20, 16, 12, 11,  5, ...
      89, 90, 95, 93, 82, 80, 65, 17, 18, 13, 14, 15, 10,  9, ...
     104,103,111,107,108, 84, 78, 77, 79, 58, 47, 48, 46, 45, 49, ...
     102,101,110, 94, 73, 81, 83, 36, 57, 42, 40, 38, 56, 51, ...
     116,115,112, 96, 75, 76, 33, 59, 43, 44, 39, 55, 50, ...
     113,114,109,  4, 74, 34, 60, 41, 37, 53, 52,...
     120,118,  3,  1, 35, 61, 63, 54, ...
     119,117,  2, 64, 62];
   
   lid2 = [105, 106]; % These two scalars may be flipped. Verify.

  % Aggregate the tile ID numbers of all tiles in the spherical arena in
  % the correct order, taking into account where exactly the lid is placed.
  switch lid
    case 'none'
      neworder = [hemisphere1,hemisphere2];
    case 'top only'
      neworder = [lid1,hemisphere1,lid2,hemisphere2];
    case 'bottom only'
      neworder = [hemisphere1,lid1,hemisphere2,lid2];
  end
  
  % The following line ensures that unassigned tile numbers are padded with
  % NaNs. If you remove the line, they will be padded with zeros instead -
  % or skipped entirely if there are no higher, actually assigned numbers.
  % To safely pad the array with zeros, replace NaN(...) with zeros(...).
  numframe = size(original,3);
  reordered = NaN(8,8,numframe,240);
  
  % Rearrange the fourth dimension of the arrays containing LED positions.
  % Remember that dimension 1 are the actual coordinates (e.g., azimuth and
  % elevation), dimensions 2 and 3 cluster the 8*8 individual LED on each
  % tile, and dimension 4 goes through all tiles in the setup. The old
  % order went through all tiles top-to-bottom, meridian-to-side; the new
  % order goes through all tiles from the tile with ID no. 1 to the tile
  % with ID no. 236 (or whatever else the maximum is).
  numtile = numel(neworder);
%   numtile = size(original,4);
  oldorder = 1:numtile;
  reordered(:,:,:,neworder) = original(:,:,:,oldorder);
  
  % Display how many tile IDs were found, and how many more could be used.
%   unassigned = numel(find(isnan(reordered)))/(64*size(original,1));
  unassigned = 240 - numtile;
  display(['-- Out of 240 supported tiles, ',num2str(numtile), ...
           ' tile IDs were assigned; ',num2str(unassigned), ...
           ' were left unassigned. --'])

  

  % SECOND, REARRANGE TILES INTO A VIRTUAL, RECTANGULAR PATTERN
  %
  % As the existing code driving our stimulus arena expects LED tiles to be
  % arranged on the surface of a cylinder (or on a flat rectangular
  % surface), we need to rearrange our distribution of LED tiles into such 
  % a rectangular array. The resulting position of certain tiles may seem
  % nonsensical (a tile from the top ending up in the centre of the
  % rectangle, its neighbour up in the bottom right corner etc.), but this
  % new rearrangement is only virtual and has no deeper meaning besides
  % getting the code to work properly. Don't worry too much about it.

  % Arrange increasing tile ID top-to-bottom, then left-to-right, in
  % vertical columns of 8. The total number of columns is 15 for 120 tiles,
  % 30 for 240 tiles.

  rearranged = NaN(64,240,numframe);
  
  % This is a hack. Clean up in the near future.
  reordered = permute(reordered,[2 1 3 4]);
  
  % This is not a hack. The following loop must always run to 240.
%   for tileID = 1:numtile
  for tileID = 1:240
    
    xshift = 8*floor((tileID-1)/8);
    yshift = 8*mod(tileID-1,8);
%     display([tileID xshift/8 yshift/8])
   
%     if tileID <= 5
%         % THIS IS A HACK FOR DEBUGGING. REMOVE.
%       rearranged((1:8)+56-yshift,(1:8)+xshift,:)  ...
%         = repmat(zebra,[1 1 numframe 1]);
%     else
      rearranged((1:8)+56-yshift,(1:8)+xshift,:) = reordered(:,:,:,tileID);
%     end
    
  end
  
  
end

function [sphereradius, ...
          ribbonradius, ribbonangle, ribangle] = sphereShape(lid)

  numribbon = 11;
  tilewidth = 21; % Here, used only to compute elevation of ribs/ribbons.
  ribwidth = 2.1; % Here, used only to compute elevation of ribs/ribbons.
  % Elsewhere, tilewidth is 20 (which is correct). Suggestion: Set ribwidth
  % to 3.1 instead, and use tilewidth = 20 throughout.
  
  % To minimize the gap near the poles of the sphere, determine the radius
  % for which the sum of the elevation angles covered by LED tiles, covered
  % by structural ribs, and covered by the desired holes at the top and
  % bottom equals approaches 180 degrees. The following equation is a cost
  % function penalising any deviation from this optimal radius.
  eqn = @(radius) abs((numribbon)   * 2*asind((tilewidth/2)/radius) + ...
                      (numribbon+1) * 2*asind((ribwidth/2)/radius) + ...
                      1             * 2*asind(60.1/(2*radius)) - 180);
  
  % Numerically solve the equation to find the optimal sphere radius.
  sphereradius = fminsearch(eqn,50);
  
  % Next, Julian decided to deviate from the optimum for practical reasons.
  % The original code did NOT work for any stretchfactors other than 1.05,
  % but this problem has been fixed since. All factors > 1.05 should work.
  stretchfactor = 1.05;
%   stretchfactor = 1;
  sphereradius = sphereradius * stretchfactor;
  
  % Here, "sphereradius" refers to the radius of the sphere on which the
  % centres of all LED tiles are located. The following "outerradius" is
  % the radius of the sphere on which the inner edges of tiles are located.
  outersphereradius = sqrt(sphereradius^2+(tilewidth/2)^2);
  
  % (Re-)Compute the elevation angles covered by a tile, and by a rib.
  ribbonangle = 2*asind((tilewidth/2)/sphereradius);
  ribangle    = 2*asind((ribwidth/2)/sphereradius);
  
  % Compute the radii of perfectly horizontal planes containing one rib
  % each. This computation assumes an odd number of ribbons, i.e. the
  % presence of an equatorial ribbon, rather than an equatorial rib.
  % Because of symmetry, only the unique radii are computed, then
  % replicated once for their mirror-symmetric counterpart.
  
  % First, compute the number of unique ribs.
  numuniquerib = (numribbon-1)/2+1;
  
  % Second, compute the unique circular radii.
  ribbonradius = NaN(numuniquerib,1);
%   Position_connector = NaN((numribbon-1)/2+1,1);

  for k = 1:((numribbon-1)/2+1);
    ribbonradius(k) = outersphereradius * cosd((k-1/2)*ribbonangle + ...
                                                    (k-1)*ribangle);
%     Position_connector(i) = outerradius * sind(ribbonangle/2 + ...
%                             (i-1) * ribbonangle + (i-1) * ribangle);
  end

  % Third, replicate the mirror-symmetric copies.
  ribbonradius = [fliplr(ribbonradius(2:numuniquerib)'), ...
                       ribbonradius'];
                     
	% Fourth, if a lid is present, add an additional circular radius.
  if ~strcmp(lid,'none')
    
    lidradius = outersphereradius * ...
                cosd((numuniquerib+.5)*ribbonangle + numuniquerib*ribangle);

    % Depending on where the lid is located, place its radius
    % at the top or at the bottom of the list of circular radii.
    switch lid
      case 'top only'
        ribbonradius = [lidradius ribbonradius];
      case 'bottom only'
        ribbonradius = [ribbonradius lidradius];
      otherwise
        error(['-- How many lids are covering the poles of the ', ...
               'sphere? Three scenarios are possible: ''none'', ', ...
               '''top only'', and ''bottom only''. Please select ', ...
               'one of them.--']);
    end
    
  end
  
  % The following line looks useless.
%   Value_for_triangle = sind(Theta_1/2)*6;

end

function tilepos = sphereTilePosition(ribbonradius, ribbonangle, ...
                                      ribangle, numtileperribbon, lid)

  ribbonbeta = (ribbonangle+ribangle) * (5:-1:-5);

  switch lid
    case 'none'
      % Relax. Do nothing.
    case 'top'
      ribbonbeta = [ribbonangle*6 + ribangle*6, ribbonbeta];
    case 'bottom'
      ribbonbeta = [ribbonbeta, -(ribbonangle*6 + ribangle*6)]; 
    otherwise
      error(['-- How many lids are covering the poles of the sphere? ', ...
             'Three scenarios are possible: ''none'', ''top only'', ', ...
             'and ''bottom only''. Please select one of them. --']);   
  end


  tilepos = NaN(sum(numtileperribbon),5);
  tilepos(:,1) = 1; % Orientation 1 = Upside down, 0 = Normal
  counter = 0;

  for a = 1:length(ribbonradius)%go through all rows
    for b = 1:numtileperribbon(a)%go through all elements of each row
      
      % Compute the azimuth of each tile centre,
      % taking into account the 6mm wide meridian "keel" or "spine",
      % and the 20mm width of each tile.
      keelwidth = 6;
      tilewidth = 20;

      % Compute the azimuth covered by the keel, and by each tile.
      % These numbers are exact for the keel, as well as for tiles along
      % the equatorial ribbon. For all other ribbons, they are ONLY
      % APPROXIMATE, because these tiles are perpendicular to the
      % sphere, but not perpendicular to the circle formed by the ribbon
      % (i.e., they are not perfectly vertical). Their true azimuth spread
      % would be slightly larger.
      keelangle = 2*asind((keelwidth/2)/ribbonradius(a));
      tileangle = 2*asind((tilewidth/2)/ribbonradius(a));
  
      % Compute the azimuth of the centre of this tile,
      % considering the azimuth offset and the tiles already placed.
      thistilealpha = keelangle + (b-1/2)*tileangle;
 
      % Counter to save in the desired order.
      counter = counter + 1; 
      tilepos(counter,3) = thistilealpha;
      tilepos(counter,4) = ribbonbeta(a);
      tilepos(counter,5) = ribbonradius(a);
      
    end

  end

end

function [ledposcartesian, ...
          ledposgeographic] = sphereLedPosition(tilepos,sphereradius)

  ledseparation = 2.48;
  numtile = size(tilepos,1);
  tileisflipped = tilepos(:,1);
  
  % Pre-allocate variable size to speed up computation.
  ledpos           = NaN(length(tilepos)*64,5);  
  ledposcartesian  = NaN(3,8,8,numtile);
  ledposgeographic = NaN(2,8,8,numtile);
  
  for tile = 1:numtile

    % Read out the position of the tile centre in geographic coordinates.
    tilealpha  = tilepos(tile,3);
    tilebeta   = tilepos(tile,4);
    tileradius = sphereradius;
    
    % Convert the position of the tile centre into cartesian coordinates.
    tilecentre = tileradius * [cosd(tilebeta)*cosd(tilealpha); ...
                               cosd(tilebeta)*sind(tilealpha); ...
                               sind(tilebeta)];

    % Compute the "unit" vector in the beta direction (still cartesian).
    betaunitvector = [-sind(tilebeta)*cosd(tilealpha); ...
                      -sind(tilebeta)*sind(tilealpha); ...
                       cosd(tilebeta)];
                     
    % Normalise the vector to make sure it is a unit vector.
    betaunitvector = betaunitvector/norm(betaunitvector);

    % Compute the "unit" vector in the alpha direction  (still cartesian).
    alphaunitvector = [-cosd(tilebeta)*sind(tilealpha); ...
                        cosd(tilebeta)*cosd(tilealpha); ...
                        0];
                      
    % Normalise the vector to make sure it is a unit vector.
    alphaunitvector = alphaunitvector/norm(alphaunitvector);

    
    % Now, place 64 individual LEDs around the tile centre.
    for row = 1:8    % go through columns (!)
      for col = 1:8  % go through rows (!)

        % Compute LED position in cartesian coordinates.
        b = ledseparation * (row - 4.5);
        a = ledseparation * (col - 4.5);
        rled = tilecentre + b*betaunitvector + a*alphaunitvector;
        
        % Save them for later.
        ledposcartesian(:,col,row,tile) = rled;

        % Convert LED position from cartesian to geographic coordinates.
        beta  = asind(rled(3)/norm(rled));
        alpha = atand(rled(2)/rled(1));

        % Constrain geographic coordinate values to standard range.
        % beta = mod(beta,180)-90; % Not needed.
        alpha = mod(alpha,180);
        
        % Save them for later.
        ledposgeographic(:,col,row,tile) = [beta,alpha];

        % Remember which tile this LED is on (i.e., its ID number).
        id = (tile-1)*64 + (row-1)*8 + col;
        
        ledpos(id,2) = tilepos(tile,2); % Remember the tile ID number.

      end
    end
  end
  
  allisflipped = floor(sum(tileisflipped)/numel(tileisflipped));
  
  if allisflipped
    
    % Were all tiles accidentally flipped upside down during construction?
    % If so, flip them back the way they belong.
%     ledposcartesian = flipdim(flipdim(ledposcartesian, 2), 3);
%     ledposgeographic = flipdim(flipdim(ledposgeographic, 2), 3);
%% REVISIT THIS!

  else
    
    % Were individual tile flipped upside down during construction?
    % Or, more precisely, rotated 180 degrees around its centre point?
    % If so, flip the LED coordinates back the way they belong.
    
    for tile = 1:numtile
      if tileisflipped(tile)

        ledposcartesian(:,:,:,tile) = ...
          flipdim(flipdim(ledposcartesian(:,:,:,tile), 2), 3);

        ledposgeographic(:,:,:,tile) = ...
          flipdim(flipdim(ledposgeographic(:,:,:,tile), 2), 3);

      end
    end
    
  end
  
end