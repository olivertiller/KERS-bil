// Diagnoseversjon av utsnurring_logger: logger begge flankene fra US1881-latchen,
// over USB og trådløst over BLE.
//
// USB (når en PC er koblet til):
//   <t_s>,S   utgangen gikk lav  (sensoren så sydpolen)
//   <t_s>,N   utgangen gikk høy  (sensoren så nordpolen)
//   # t=<s> nivaa=<L/H> S=<antall> N=<antall> tapt=<antall>   statuslinje hvert sekund
//
// BLE: enheten heter "KERS-utsnurring" og har én tjeneste med to karakteristikker.
//   FLANKER  notify, binært: n poster à 5 byte, uint32 micros() (little endian) + 'S'/'N'
//   STATUS   notify, samme statuslinje som over, hvert sekund
// Flankene ligger i en ringbuffer på 8192 poster til de er sendt. Faller forbindelsen
// kort ut, sendes det som kom imens når den er tilbake. Ved ny tilkobling sendes bare
// flanker fra de siste 30 sekundene, så gamle målinger ikke dukker opp i en ny serie.
//
// Riktig oppførsel er S og N annenhver gang. To like etter hverandre betyr at sensoren
// gikk glipp av den andre magneten. Kobling: VDD -> 5V, GND -> GND, OUT -> D9 (utsnurring_logger bruker D10).
// Les med  python3 Code/utsnurring_diagnose/les_diagnose.py [--ble]

#include <BLEDevice.h>
#include <BLEServer.h>
#include <BLE2902.h>

const uint8_t HALL_PIN = D9;
const bool LED_AKTIV_LAV = true;

const char *NAVN = "KERS-utsnurring";
const char *TJENESTE_UUID = "eedccbc9-a16a-4a60-bd93-a2a3324f97dc";
const char *FLANKER_UUID = "4a4bf346-512e-49b7-97a3-f5c3cd2b8a71";
const char *STATUS_UUID = "c0ebe3b8-889d-414e-87a1-23cde9e22c61";

const uint32_t BUF = 8192;              // poster, 8192 * 5 byte = 40 kB
const uint32_t MAKS_ALDER_US = 30000000; // eldre flanker sendes ikke ved ny tilkobling
const uint8_t POST = 5;
const uint8_t MAKS_POSTER = 48;          // per notify, 240 byte

volatile uint32_t tider[BUF];
volatile uint8_t nivaaer[BUF];
volatile uint32_t skrevet = 0;  // totalt antall flanker, indeks = skrevet % BUF
volatile uint32_t antall_s = 0, antall_n = 0;

uint32_t les_serie = 0, les_ble = 0, tapt = 0;
unsigned long sist_status = 0, sist_sendt = 0;

BLEServer *server;
BLECharacteristic *flanker, *status;
volatile bool tilkoblet = false, ny_tilkobling = false;

void IRAM_ATTR flanke() {
  uint8_t n = digitalRead(HALL_PIN);
  uint32_t i = skrevet % BUF;
  tider[i] = micros();
  nivaaer[i] = n;
  if (n) antall_n++; else antall_s++;
  skrevet++;
}

class Tilkobling : public BLEServerCallbacks {
  void onConnect(BLEServer *s) override {
    tilkoblet = true;
    ny_tilkobling = true;
  }
  void onDisconnect(BLEServer *s) override {
    tilkoblet = false;
    BLEDevice::startAdvertising();
  }
};

// Henter post nr. k. Returnerer false hvis den er overskrevet.
bool hent(uint32_t k, uint32_t &t, uint8_t &n) {
  noInterrupts();
  bool ok = skrevet - k <= BUF;
  t = tider[k % BUF];
  n = nivaaer[k % BUF];
  interrupts();
  return ok;
}

void send_ble() {
  uint32_t s = skrevet;
  if (s - les_ble > BUF) {  // bufferen gikk rundt før vi fikk sendt
    tapt += s - les_ble - BUF;
    les_ble = s - BUF;
  }
  if (ny_tilkobling) {  // hopp over flanker fra før denne tilkoblingen
    ny_tilkobling = false;
    uint32_t t; uint8_t n;
    while (les_ble != s && hent(les_ble, t, n) && micros() - t > MAKS_ALDER_US) les_ble++;
  }
  if (les_ble == s) return;

  // Send når det er nok til en full pakke, eller når det har gått 50 ms.
  uint16_t mtu = server->getPeerMTU(server->getConnId());
  uint8_t maks = constrain((mtu > 3 ? mtu - 3 : 20) / POST, 1, MAKS_POSTER);
  if (s - les_ble < maks && millis() - sist_sendt < 50) return;

  uint8_t pakke[MAKS_POSTER * POST];
  uint8_t antall = 0;
  while (les_ble != s && antall < maks) {
    uint32_t t; uint8_t n;
    if (!hent(les_ble, t, n)) break;
    uint8_t *p = pakke + antall * POST;
    p[0] = t; p[1] = t >> 8; p[2] = t >> 16; p[3] = t >> 24;
    p[4] = n ? 'N' : 'S';
    antall++;
    les_ble++;
  }
  flanker->setValue(pakke, antall * POST);
  flanker->notify();
  sist_sendt = millis();
}

void setup() {
  Serial.begin(115200);
  Serial.setTxTimeoutMs(0);  // ikke vent på USB når ingen leser porten, da henger BLE etter
  pinMode(HALL_PIN, INPUT_PULLUP);
  pinMode(LED_BUILTIN, OUTPUT);

  BLEDevice::init(NAVN);
  BLEDevice::setMTU(247);
  server = BLEDevice::createServer();
  server->setCallbacks(new Tilkobling());
  BLEService *tjeneste = server->createService(TJENESTE_UUID);
  flanker = tjeneste->createCharacteristic(FLANKER_UUID, BLECharacteristic::PROPERTY_NOTIFY);
  flanker->addDescriptor(new BLE2902());
  status = tjeneste->createCharacteristic(STATUS_UUID, BLECharacteristic::PROPERTY_NOTIFY | BLECharacteristic::PROPERTY_READ);
  status->addDescriptor(new BLE2902());
  tjeneste->start();
  BLEAdvertising *reklame = BLEDevice::getAdvertising();
  reklame->addServiceUUID(TJENESTE_UUID);
  reklame->setScanResponse(true);
  BLEDevice::startAdvertising();

  attachInterrupt(digitalPinToInterrupt(HALL_PIN), flanke, CHANGE);
  Serial.println("t_s,pol");
}

void loop() {
  bool sydpol = !digitalRead(HALL_PIN);
  digitalWrite(LED_BUILTIN, sydpol != LED_AKTIV_LAV);

  // USB: skriv alt som er nytt. Uten PC tilkoblet hoppes det over.
  uint32_t s = skrevet;
  if (s - les_serie > BUF) les_serie = s - BUF;
  while (les_serie != s) {
    uint32_t t; uint8_t n;
    bool ok = hent(les_serie, t, n);
    les_serie++;
    if (!ok || !Serial) continue;
    char linje[32];
    sprintf(linje, "%lu.%06lu,%c", t / 1000000UL, t % 1000000UL, n ? 'N' : 'S');
    Serial.println(linje);
  }

  if (tilkoblet) send_ble();

  unsigned long na = millis();
  if (na - sist_status >= 1000) {
    sist_status = na;
    char linje[80];
    sprintf(linje, "# t=%lu nivaa=%c S=%lu N=%lu tapt=%lu", na / 1000,
            digitalRead(HALL_PIN) ? 'H' : 'L', antall_s, antall_n, tapt);
    if (Serial) Serial.println(linje);
    if (tilkoblet) {
      status->setValue(linje);
      status->notify();
    }
  }
}
