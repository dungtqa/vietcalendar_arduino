#ifndef VietCalendar_h
#define VietCalendar_h

#include "Arduino.h"

#define PI 3.1415926535897932384626433832795

class VietCalendar {
  public:
    long jdFromDate(int dd, int mm, int yy);
    int* jdToDate(long jd);
    double SunLongitude(double jdn);
    double SunLongitudeAA98(double jdn);
    double NewMoon(long k);
    double NewMoonAA98(long k);
    long INT(double d);
    double getSunLongitude(long dayNumber, double timeZone);
    long getNewMoonDay(long k, double timeZone);
    long getLunarMonth11(int yy, double timeZone);
    int getLeapMonthOffset(int a11, double timeZone);
    int* convertSolar2Lunar(int dd, int mm, int yy, double timeZone);
    int* convertLunar2Solar(int lunarDay, int lunarMonth, int lunarYear, int lunarLeap, double timeZone);
};

#endif
