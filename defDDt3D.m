function [D,Dt] = defDDt3D
        % defines finite difference operator D
        % and its transpose operator
        
        D = @(U) ForwardD(U);
        Dt = @(DU) Dive(DU);
        
         function DU = ForwardD(U)
            DU = cell(3,1);
            % Forward finite difference operator
            DU{1} = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
            DU{2} = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
            temp = zeros(size(U));
            temp(:,:,1:end-1) = diff(U,1,3);
            temp(:,:,end) = U(:,:,1) - U(:,:,end);
            DU{3} = temp;
            clear temp;
        end
        
        function DtXY = Dive(DU)
            X = DU{1};
            Y = DU{2};
            t = DU{3};
            % Transpose of the forward finite difference operator
            DtXY = [X(:,end,:) - X(:, 1,:), -diff(X,1,2)];
            DtXY = DtXY + [Y(end,:,:) - Y(1, :,:); -diff(Y,1,1)];
             temp = zeros(size(t));
             temp(:,:,2:end) = -diff(t,1,3);
             temp(:,:,1) = t(:,:,end) - t(:,:,1);
             DtXY = DtXY + temp;
        end
        
    end