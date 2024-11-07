public class Main {
    public static void main(String[] args) {
        ShootingMethod shootingMethod = new ShootingMethod(1.0, 3.0, 25, 0.001, 10000, -25.0);
        shootingMethod.bisectionMethod();
    }
}
