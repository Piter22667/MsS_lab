package test3;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import java.awt.BasicStroke;
import java.awt.Color;
import java.io.File;
import java.io.IOException;

public class SystemResponsePlotter {
    private static final int WIDTH = 800;
    private static final int HEIGHT = 600;

    public static void createPlot(double[][] measurements, double[][] calculated, String filename) {
        XYSeriesCollection dataset = new XYSeriesCollection();


        Color[] colors = {
                new Color(0, 0, 255),    // Синій для y1
                new Color(255, 0, 0),    // Червоний для y2
                new Color(0, 100, 0)     // Зелений для y3
        };

        for (int i = 0; i < 3; i++) {
            XYSeries measuredSeries = new XYSeries("Виміряні y" + (i+1));
            XYSeries modelSeries = new XYSeries("Модель y" + (i+1));

            // Додаємо точки
            for (int j = 0; j < measurements[0].length; j++) {
                double t = j * 0.2; // Час (DT = 0.2)
                measuredSeries.add(t, measurements[i*2][j]);
                modelSeries.add(t, calculated[i*2][j]);
            }

            dataset.addSeries(measuredSeries);
            dataset.addSeries(modelSeries);
        }

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Відгук системи",
                "Час",
                "Амплітуда",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYPlot plot = chart.getXYPlot();
        plot.setDrawingSupplier(null);  // Скидаємо налаштування за

        for (int i = 0; i < 3; i++) {
            int measuredIndex = i * 2;
            int modelIndex = i * 2 + 1;

            plot.getRenderer().setSeriesPaint(measuredIndex, colors[i]);
            plot.getRenderer().setSeriesStroke(measuredIndex,
                    new BasicStroke(2.0f));

            plot.getRenderer().setSeriesPaint(modelIndex, colors[i]);
            float[] dash = {5.0f, 5.0f};
            plot.getRenderer().setSeriesStroke(modelIndex,
                    new BasicStroke(
                            2.0f,                    // Товщина лінії
                            BasicStroke.CAP_SQUARE,
                            BasicStroke.JOIN_MITER,
                            10.0f,                   // Масштаб штрихів
                            dash,                    // Патерн штрихів
                            0.0f                     // Фаза пунктиру
                    )
            );
        }

        plot.setDomainGridlinesVisible(true);
        plot.setRangeGridlinesVisible(true);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlineStroke(new BasicStroke(0.5f));
        plot.setRangeGridlineStroke(new BasicStroke(0.5f));

        // Встановлюємо діапазони осей
        plot.getDomainAxis().setRange(0, 52);
        plot.getRangeAxis().setRange(-1.3, 1.2);

        try {
            ChartUtils.saveChartAsPNG(new File(filename), chart, WIDTH, HEIGHT);
        } catch (IOException e) {
            System.err.println("Помилка при збереженні графіку: " + e.getMessage());
        }
    }
}