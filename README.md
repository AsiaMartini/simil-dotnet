# simil-dotnet (NOT YET WORKING)

.NET porting of https://github.com/gabriel-de-luca/simil

## Usage

```
float[,] local_ctrl_p = new float[,] {
    {0, 0, 0},
    {0, 12.24f, 0},
    {32.42f, 6.79f, 0}
};

GeodeticToEcef(43.88072176381581f, 11.072238294111774f, 53, out double x1, out double y1, out double z1);
GeodeticToEcef(43.88072176381581f, 11.072238294111774f, 53, out double x2, out double y2, out double z2);
GeodeticToEcef(43.88072176381581f, 11.072238294111774f, 53, out double x3, out double y3, out double z3);

float[,] geocent_ctrl_p = new float[,] {
    {float.Parse(x1.ToString()), float.Parse(y1.ToString()), float.Parse(z1.ToString())},
    {float.Parse(x2.ToString()), float.Parse(y2.ToString()), float.Parse(z2.ToString())},
    {float.Parse(x3.ToString()), float.Parse(y3.ToString()), float.Parse(z3.ToString())},
};

float m;
NDarray r, t;
(m, r, t) = Simil.process(local_ctrl_p, geocent_ctrl_p);
```
