function [half_axis_big, half_axis_small] = get_axis(image_binor)
    %����������� ������� ����� - ����� � ������� ���� ������
    % �� ���� �������� �������
    stats = regionprops('table',image_binor,'MajorAxisLength','MinorAxisLength');
    half_axis_big = stats.MajorAxisLength/2;
    half_axis_small = stats.MinorAxisLength/2;
end  