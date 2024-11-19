package test3;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import java.awt.BasicStroke;
import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

public class ConvergencePlotter {
    private static final int WIDTH = 800;
    private static final int HEIGHT = 600;

    public static void createPlot(List<Double> qualityHistory, String filename) {
        XYSeries series = new XYSeries("Збіжність");

        for (int i = 0; i < qualityHistory.size(); i++) {
            series.add(i, qualityHistory.get(i));
        }

        // Створюємо набір даних
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);

        // Створюємо графік
        JFreeChart chart = ChartFactory.createXYLineChart(
                "Збіжність",           // Заголовок
                "Ітерація",           // Підпис осі X
                "Показник якості",    // Підпис осі Y
                dataset              // Дані
        );

        chart.getXYPlot().getRenderer().setSeriesPaint(0, Color.BLUE);
        chart.getXYPlot().getRenderer().setSeriesStroke(0, new BasicStroke(2.0f));

        chart.getXYPlot().getDomainAxis().setRange(-0.1, qualityHistory.size() - 0.9);


        try {
            // Зберігаємо графік у файл
            ChartUtils.saveChartAsPNG(
                    new File("convergence.png"),
                    chart,
                    WIDTH,
                    HEIGHT
            );
        } catch (IOException e) {
            System.err.println("Помилка при збереженні графіку: " + e.getMessage());
        }
    }
}