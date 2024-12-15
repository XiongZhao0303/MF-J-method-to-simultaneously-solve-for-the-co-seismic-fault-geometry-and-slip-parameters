function L = calculateL(Faults, f)
% calculateL - Computes the Laplacian matrix for an arbitrary number of faults.
%
% Inputs:
%   Faults - A cell array where each cell contains the fault matrix for one fault.
%   f - A cell array where each cell contains the [ast, adi] for the corresponding fault.
%
% Output:
%   L - The combined Laplacian matrix for all faults.

% Initialize the combined Laplacian matrix

L = [];

% Iterate over each fault
for faultIdx = 1:length(Faults)
    Fault = Faults{faultIdx};
    ast = f{faultIdx}(1);
    adi = f{faultIdx}(2);
    
    xc = Fault(:, 1);
    len = Fault(:, end-2);
    width = Fault(:, end-1);
    
    NW = adi - 1;
    NL = ast - 1;
    T = [];
    k = 1; % Index for T matrix
    
    % Compute the Laplacian matrix for the current fault
    for j = 1:NW
        for i = 1:NL
            for m2 = 1:2
                index1 = (j-1)*NL + i;
                index2 = (j-1)*NL + i - 1;
                index3 = (j-1)*NL + i + 1;
                index4 = (j-2)*NL + i;
                index5 = j*NL + i;
                dx = len(1);
                dy = width(1);
                
                % Fill the Laplacian matrix entries
                if (index1 >= 1 && index1 <= length(xc))
                    T(k, 2*(index1-1)+m2) = -2*(dx.^-2 + dy.^-2);
                end
                if (index2 >= 1 && index2 <= length(xc))
                    T(k, 2*(index2-1)+m2) = dx.^-2;
                end
                if (index3 >= 1 && index3 <= length(xc))
                    T(k, 2*(index3-1)+m2) = dx.^-2;
                end
                if (index4 >= 1 && index4 <= length(xc))
                    T(k, 2*(index4-1)+m2) = dy.^-2;
                end
                if (index5 >= 1 && index5 <= length(xc))
                    T(k, 2*(index5-1)+m2) = dy.^-2;
                end
                k = k + 1;
            end
        end
    end
    
    % Append the current fault's Laplacian matrix to the combined matrix
    L = blkdiag(L, T);
end
end
