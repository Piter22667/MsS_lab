package test3;

public class Result {
    public final ParameterIdentification.SystemParameters parameters;
    public final double quality;
    public final int iterations;

    public Result(ParameterIdentification.SystemParameters parameters, double quality, int iterations) {
        this.parameters = parameters;
        this.quality = quality;
        this.iterations = iterations;
    }
}