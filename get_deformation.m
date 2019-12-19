function [VAR] = get_deformation(hor_ax_small, vert_ax_big, VAR)
    l = vert_ax_big*VAR.pixel_size; % ������� �� ����� � �������� � ��
    d = hor_ax_small*VAR.pixel_size;
    
%     VAR.R = (l+d)/2; % ����� ����������� ������� �� ������ ���������� 
%     Volume = d^2*l; % � ������� ��������

    VAR.R = (d.^2.*l).^(1/3); % ����� ����������� ������� ����� - ��� ������� 
    % ���� � ������� ������ ������ �������. �������� ��� � ��������������� 
    % ����������, ��� � � ������ ���������� ��������� ����������
    
    VAR.D = (l-d)./(l+d)*100; %������ ����������
    VAR.time = linspace(VAR.frame_period,length(VAR.D)*VAR.frame_period,length(VAR.D)); 
    %����� ����������
end
