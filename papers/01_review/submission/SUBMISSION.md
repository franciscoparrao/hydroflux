# Submission a EMS — lista de verificación

**Generado**: 2026-08-18 · **Commit pinneado en el paper**: `ca0062e`
**Estado**: `preflight.sh` PASSED en los 15 chequeos.

Sistema: Editorial Manager de Elsevier, vía la página del journal.
La subida es manual — esta lista es para seguirla en pantalla.

---

## Archivos de este directorio

| Archivo | Qué es | Item en Editorial Manager |
|---|---|---|
| `manuscript.pdf` | 53 pp, elsarticle, figuras y tablas embebidas | Manuscript |
| `highlights.txt` | 5 bullets, máx 84 de 85 caracteres | Highlights |
| `graphical_abstract.pdf` | 3070 × 1181 px (mín. 1328 × 531) | Graphical Abstract |
| `cover_letter.pdf` | 3 pp | Cover Letter |

Los cuatro están regenerados desde el manuscrito actual. Ninguno arrastra
texto de versiones anteriores — verificado con grep sobre el texto
extraído de cada PDF.

---

## Datos que el sistema va a pedir

**Título**
hydroflux: A Well-Balanced, Differentiable-by-Design 2D Shallow-Water
Solver in Rust, Verified Against Analytical and Community Benchmarks and
Applied to a Semiarid Andean Reach

**Tipo de artículo**: Research Article

**Autores, en este orden** (el orden debe coincidir con el del manuscrito):

| # | Autor | Afiliación | ORCID |
|---|---|---|---|
| 1 | Francisco Parra (autor de correspondencia) | Depto. Ing. Informática, USACH, Santiago, Chile | 0009-0006-0435-1854 |
| 2 | Verónica Gil-Costa | UNSL y CONICET, San Luis, Argentina | 0000-0003-4637-9725 |
| 3 | Carolina Bonacic | USACH, Santiago, Chile | 0000-0002-8076-6537 |
| 4 | Mauricio Marín | USACH, Santiago, Chile | 0000-0003-0662-7149 |

> **Ojo con el ORCID del autor 1**: usar `0009-0006-0435-1854`, que es el
> del paper de SurtGIS publicado. El otro registro
> (`0009-0008-4961-304X`) es un duplicado; conviene fusionarlo en ORCID
> en algún momento, pero no antes de someter.

**Keywords** (1–7): shallow water equations, finite volume,
well-balanced, differentiable physics, Rust, flood modelling

**Grant/contract number**: `062619MC_POSTDOC`

Es el código del proyecto DICYT. El numeral 10.7 de las bases exige que
los agradecimientos digan "Agradecimientos al Proyecto ___, Código ___,
Vicerrectoría de Investigación, Innovación y Creación", y que la
afiliación del postdoctorante sea Universidad de Santiago de Chile en
exclusiva. Ambos requisitos están cumplidos en el manuscrito.

**Data statement**: el manuscrito trae sección *Data availability*, que
cumple la Option C de Elsevier con depósito real: el archivo citable
está publicado en Zenodo bajo **doi:10.5281/zenodo.22014346**
(v0.1.0, 4.4 MB), citado y enlazado en el artículo. Contiene solver,
tests, benchmarks, insumos de la aplicación y scripts de figuras;
excluye las fuentes del manuscrito.

**Declaración de conflicto de interés**: en el manuscrito, sin conflictos.

**Declaración de uso de IA**: en el manuscrito y en el cover letter.

**Revisores sugeridos**: en el cover letter (Caviedes-Voullième,
Roberts, García-Navarro, Shen, Escauriaza).

---

## Antes de apretar submit

- [ ] Correr `bash papers/01_review/preflight.sh` — debe decir PASSED.
      Ha cazado tres veces que el pin quedó atrasado tras un commit.
- [ ] Confirmar que `origin/main` está al día: el paper apunta a un
      commit público y si no está pusheado, la instrucción de
      reproducción falla.
- [ ] Revisar el cover letter una vez con ojo propio. Es lo único que
      el editor lee antes de decidir si manda a revisión.

## Pendientes que NO bloquean

- Agregar el DOI del artículo a `related_identifiers` del registro de
  Zenodo cuando exista. Zenodo permite editar metadata de un registro
  publicado sin cambiar el DOI, así que es un trámite posterior.
- Dirección postal completa en las afiliaciones. Hoy llevan
  organización, ciudad y país, que es el formato estándar de
  elsarticle; las guidelines piden calle y código postal.
- El "más de un orden de magnitud" de §3.4 reemplazó a un "45×" que ya
  no es verificable, porque la variante contra la que comparaba no
  existe en el código.
