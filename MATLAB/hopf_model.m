classdef hopf_model < handle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                               LICENSE                             %%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Copyright 2019 Mario Senden
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                             DESCRIPTION                           %%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Z = hopf_model(); creates an instance of the Hopf model using  
% standard parameter values see: Senden, M., Reuter, N., 
% van den Heuvel, M. P., Goebel, R., & Deco, G. (2017). 
% Cortical rich club regions can organize state-dependent 
% functional network formation by engaging in oscillatory behavior. 
% NeuroImage, 146, 561–574. 
% https://doi.org/10.1016/j.neuroimage.2016.10.044
% 
% optional inputs are
% 
% - connectivity: structural connectivity matrix
% - omega       : intrinsic frequency of each region (in radians)
% - alpha       : bifurcation parameter
% - g_coupling  : global coupling strength
% - sigma       : noise scaling factor
% - tau         : time constant (in seconds)
% - dt          : integration time step (in seconds)
% 
% 
% this class has the following functions
% 
% - signal = hopf_model.simulate(time,sampling_rate)
% - [X,Y] = hopf_model.get_activity()
% - hopf_model.relax(time)
% - hopf_model.set_connectivity(connectivity)
% - hopf_model.set_omega(omega)
% - hopf_model.set_alpha(alpha)
% - hopf_model.set_coupling(g_coupling)
% - hopf_model.set_tau(tau)
% - hopf_model.set_dt(dt)

    properties (Access = private)
        % functions
        
        % model parameters
        n_regions                       % number of regions
        C                               % connectivity matrix
        omega                           % angular velocity
        alpha                           % bifurcation parameter
        g_coupling                      % global coupling factor
        sigma                           % noise scaling factor
        tau                             % time constant
        
        % simulation parameters
        dt                              % integration time step
        dsig                            % noise integrator
        
        % outcomes
        x                               % real part
        y                               % imaginary part
        
    end
    
    methods (Access = public)
        % constructor
        function self = hopf_model(varargin)
            p = inputParser;
            addOptional(p,'connectivity',[]);
            addOptional(p,'omega',[]);
            addOptional(p,'alpha',0);
            addOptional(p,'g_coupling',0.1);
            addOptional(p,'sigma',0.02);
            addOptional(p,'tau',1);
            addOptional(p,'dt',1e-1);
            
            p.parse(varargin{:});
            
            self.C = p.Results.connectivity;
            self.n_regions = size(self.C,1);
            self.omega = p.Results.omega;
            self.alpha = p.Results.alpha;
            self.g_coupling = p.Results.g_coupling;
            self.sigma = p.Results.sigma;
            self.tau = p.Results.tau;
            self.dt = p.Results.dt;
            self.dsig = self.sigma * sqrt(self.dt / self.tau);
            
            self.x = zeros(self.n_regions,1);
            self.y = zeros(self.n_regions,1);
            
        end
        
        function [X,Y] = get_activity(self)
            X = self.x;
            Y = self.y;
        end
        
        function set_connectivity(self,C)
           self.C = C;
           self.n_regions = size(C,1);
        end
        
        function set_omega(self,omega)
            self.omega = omega;
        end
        
        function set_alpha(self,alpha)
            self.alpha = alpha;
        end
        
        function set_coupling(self,G)
            self.g_coupling = G;
        end
        
        function set_tau(self,tau)
            self.tau = tau;
        end
        
        function set_dt(self,dt)
            self.dt = dt;
        end
        
        function relax(self,time)
            t_steps = time / self.dt + 1;
            for t=1:t_steps
                self.update();
            end
        end
        
        function signal = simulate(self,time,sr)
            t_steps = time / self.dt + 1;
            t_bold = time / sr + 1;
            signal = zeros(t_bold,self.n_regions);
            idx = 1;
            for t=0:t_steps-1
                self.update();
                if ~mod(t, sr / self.dt)
                    signal(idx,:) = self.x;
                    idx = idx + 1;
                end
            end
        end
    end
    methods (Access = private)
        function update(self)
            dx = self.dt / self.tau *...
                ((self.alpha - self.x.^2 - self.y.^2) .* self.x -...
                self.omega .* self.y + self.g_coupling *...
                sum(self.C .* (meshgrid(self.x) - meshgrid(self.x)'),2)) + ...
                self.dsig * randn(self.n_regions,1);
            dy = self.dt / self.tau *...
                ((self.alpha - self.x.^2 - self.y.^2) .* self.y +...
                self.omega .* self.x + self.g_coupling *...
                sum(self.C .* (meshgrid(self.y) - meshgrid(self.y)'),2)) + ...
                self.dsig * randn(self.n_regions,1);
            
            self.x = self.x + dx;
            self.y = self.y + dy;
        end
    end
end
