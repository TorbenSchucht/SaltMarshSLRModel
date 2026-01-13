function plotVerticalChange(redistributionChange, va, sea_level_rise,h_domain, front_index, t,x)

%%
subplot(2,1,1)
plot(x, h_domain, 'k','LineWidth', 2)
ylim([-1.3, 1.3])
xlim([-1,1])
xlabel('x')
ylabel('h(x)')
yline(0, '--', 'LineWidth', 1)
xline(x(front_index), '--', 'LineWidth', 1)
text(-0.32,0.084,'Marsh edge', 'Rotation', 90, 'FontSize', 10)
text(0.4, 0.2,'MHWN', 'FontSize', 10)

subplot(2,1,2)
plot1 = plot(x, redistributionChange(t,:), 'LineWidth', 2);
hold on
xline(x(front_index), '--', 'LineWidth', 1)
plot2 = plot(x, va(t,:), 'LineWidth', 2);
yline(0, '--', 'LineWidth', 1)
ylineSea = yline(-sea_level_rise, 'r', 'LineWidth', 2, 'Color', [0.09290, 0.6940, 0.1250]);
xlabel('x')
ylabel('\Delta h')

legend([plot2, plot1, ylineSea], {'Sedimentation', 'Redistribution', 'SLR'})


end