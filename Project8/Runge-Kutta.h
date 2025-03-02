#pragma once

class DerivativeBase {
public:
  virtual float derivative(float x, float y) = 0;
};
class RungeKutta {
private:
  float k1 = 0, k2 = 0, k3 = 0, k4 = 0, h, y, x;
  DerivativeBase *derivativeObj;

public:
  RungeKutta(float h, float y, float x, DerivativeBase *derivativeObj): h((h != 0) ? h : 0.001), y(y), x(x), derivativeObj(derivativeObj) {}

  inline void calcK1() { k1 = h * derivativeObj->derivative(x, y); }

  inline void calcK2()
  {
    k2 = h * derivativeObj->derivative(x + h / 2, y + k1 / 2);
  }
  inline void calcK3() {
    k3 = h * derivativeObj->derivative(x + h / 2, y + k2 / 2);
  }
  inline void calcK4() {
      k4 = h * derivativeObj->derivative(x + h, y + k3);
  }

  void step() {
    calcK1();
    calcK2();
    calcK3();
    calcK4();
    y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    x += h;
  }
    float getY() const { return y; }
    float getX() const { return x; }

    void setY(float newY) { y = newY; }
    void setX(float newX) { x = newX; }

    void set_h(float new_h) { h = new_h; }
};
