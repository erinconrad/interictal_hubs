function out_point = convert_plot_to_figure_coords(point,curr_ax)
ax = axis(curr_ax);

Xrange=ax(2)-ax(1);
Yrange=ax(4)-ax(3); 
x = point(1);
X=(x-ax(1))/Xrange +ax(1)/Xrange;
out_point(1) = X;

if length(point) > 1
    y = point(2);
    Y=(y-ax(3))/Yrange +ax(3)/Yrange;
    out_point(2) = Y;
end