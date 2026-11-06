// Diagnoseversjon av utsnurring_logger: logger begge flankene fra US1881-latchen.
//
//   <t_s>,S   utgangen gikk lav  (sensoren så sydpolen)
//   <t_s>,N   utgangen gikk høy  (sensoren så nordpolen)
//   # t=<s> nivaa=<L/H> S=<antall> N=<antall>   statuslinje hvert sekund
//
// Riktig oppførsel er S og N annenhver gang. To like etter hverandre betyr at sensoren
// gikk glipp av den andre magneten. Kobling som utsnurring_logger (XIAO C6: OUT -> D10).
// Les med  python3 Code/utsnurring_diagnose/les_diagnose.py

#if defined(ARDUINO_ARCH_ESP32)
const uint8_t HALL_PIN = D10;
const bool LED_AKTIV_LAV = true;
#define ISR_ATTR IRAM_ATTR
#else
const uint8_t HALL_PIN = 2;
const bool LED_AKTIV_LAV = false;
#define ISR_ATTR
#endif

const uint8_t BUF = 128;

volatile unsigned long tider[BUF];
volatile uint8_t nivaaer[BUF];
volatile uint8_t skriv = 0;
volatile unsigned long antall_s = 0, antall_n = 0;
uint8_t les = 0;
unsigned long sist_status = 0;

void ISR_ATTR flanke() {
  uint8_t n = digitalRead(HALL_PIN);
  tider[skriv] = micros();
  nivaaer[skriv] = n;
  if (n) antall_n++; else antall_s++;
  skriv = (skriv + 1) % BUF;
}

void setup() {
  Serial.begin(115200);
  pinMode(HALL_PIN, INPUT_PULLUP);
  pinMode(LED_BUILTIN, OUTPUT);
  attachInterrupt(digitalPinToInterrupt(HALL_PIN), flanke, CHANGE);
  Serial.println("t_s,pol");
}

void loop() {
  bool sydpol = !digitalRead(HALL_PIN);
  digitalWrite(LED_BUILTIN, sydpol != LED_AKTIV_LAV);

  while (les != skriv) {
    noInterrupts();
    unsigned long t = tider[les];
    uint8_t n = nivaaer[les];
    interrupts();
    les = (les + 1) % BUF;

    char linje[32];
    sprintf(linje, "%lu.%06lu,%c", t / 1000000UL, t % 1000000UL, n ? 'N' : 'S');
    Serial.println(linje);
  }

  unsigned long na = millis();
  if (na - sist_status >= 1000) {
    sist_status = na;
    char linje[64];
    sprintf(linje, "# t=%lu nivaa=%c S=%lu N=%lu", na / 1000, digitalRead(HALL_PIN) ? 'H' : 'L', antall_s, antall_n);
    Serial.println(linje);
  }
}
