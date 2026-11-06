// Logger pulstider fra en US1881 Hall-latch for utsnurringstest av svinghjulet.
//
// Magneter: to stk rett overfor hverandre på svinghjulet, den ene med sydpolen og den
// andre med nordpolen mot sensoren. Det gir én puls per omdreining og et balansert hjul.
//
// Kobling (TO-92UA sett fra siden med skrift og avrundede hjørner, beina ned: 1 VDD, 2 GND, 3 OUT):
//   XIAO ESP32C6:  VDD -> 5V (VBUS)   GND -> GND   OUT -> D10 (GPIO18)
//   Arduino Uno:   VDD -> 5V          GND -> GND   OUT -> pin 2
// Utgangen er open drain og trekkes opp av den interne pull-upen. På ESP32 må en eventuell
// ekstern pull-up gå til 3V3, aldri til 5V (pinnene tåler ikke 5 V).
//
// Arduino IDE for XIAO C6: Tools > USB CDC On Boot > Enabled, ellers kommer ingen utskrift.
//
// Bruk: åpne Serial Monitor (115200 baud), spinn opp svinghjulet og la det stoppe av seg selv.
// Kopier utskriften til en .csv-fil og kjør  python3 kers_sim.py --utsnurring fil.csv
// Den innebygde LED-en lyser når sensoren ser en sydpol, nyttig for å sjekke magnetene.

#if defined(ARDUINO_ARCH_ESP32)
const uint8_t HALL_PIN = D10;
const bool LED_AKTIV_LAV = true;  // LED-en på XIAO C6 (GPIO15) lyser når pinnen er lav
#define ISR_ATTR IRAM_ATTR
#else
const uint8_t HALL_PIN = 2;
const bool LED_AKTIV_LAV = false;
#define ISR_ATTR
#endif

const byte BUF = 64;

volatile unsigned long tider[BUF];
volatile byte skriv = 0;
byte les = 0;

void ISR_ATTR puls() {
  tider[skriv] = micros();
  skriv = (skriv + 1) % BUF;
}

void setup() {
  Serial.begin(115200);
  pinMode(HALL_PIN, INPUT_PULLUP);
  pinMode(LED_BUILTIN, OUTPUT);
  attachInterrupt(digitalPinToInterrupt(HALL_PIN), puls, FALLING);
  Serial.println("t_s");
}

void loop() {
  bool sydpol = !digitalRead(HALL_PIN);  // utgangen er lav når sensoren ser en sydpol
  digitalWrite(LED_BUILTIN, sydpol != LED_AKTIV_LAV);

  // Skriv ut pulsene som har kommet siden sist. Utskriften skjer her og ikke i
  // avbruddsrutinen, slik at ingen pulser går tapt mens Serial er opptatt.
  while (les != skriv) {
    noInterrupts();
    unsigned long t = tider[les];
    interrupts();
    les = (les + 1) % BUF;

    char linje[24];
    sprintf(linje, "%lu.%06lu", t / 1000000UL, t % 1000000UL);
    Serial.println(linje);
  }
}
