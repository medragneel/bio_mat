"use strict"

const pm = document.querySelector('.pwm')
const cm = document.querySelector('.count')
const fm = document.querySelector('.freq')
const s = document.querySelector('.seqs')
const btn = document.querySelector('.btn')
const open = document.querySelector('.add')
const close = document.querySelector('.close')
const dialog = document.querySelector('dialog')

// ln function
const ln = (x) => {
    const n = 1000
    return n * (Math.pow(x, (1 / n)) - 1)
}

// log base n function
const log = (base, n) => {
    return ln(n) / ln(base)
}


// create a matrix 
const createMat = (nrows, ncols) => {
    return Array(nrows).fill().map(() => Array(ncols).fill(parseFloat(0.0)))
}

const printMat = (mat, alphabet = "ACGT") => {
    let obj = new Object()
    for (const i in alphabet) {
        obj[alphabet[i]] = mat[i]
    }
    console.table(obj)
}

const matCount = (seqs, alphabet = "ACGT") => {
    const counts = createMat(alphabet.length, seqs[0].length)
    for (const s of seqs) {
        for (let i = 0; i < seqs[0].length; i++) {
            const idx = alphabet.indexOf(s[i])
            counts[idx][i] = counts[idx][i] + 1
        }

    }

    return counts
}


const fqMat = (seqs, alphabet = "ACGT") => {
    const c = matCount(seqs)
    const f = createMat(alphabet.length, seqs[0].length)
    for (const s of seqs) {
        for (let i = 0; i < seqs[0].length; i++) {
            const idx = alphabet.indexOf(s[i])
            f[idx][i] = parseFloat((c[idx][i] / seqs.length).toFixed(2).replace(/\.?0*$/, ''));
        }

    }

    return f

}
const pwm = (seqs, alphabet = "ACGT") => {
    const f = fqMat(seqs)
    const p = createMat(alphabet.length, seqs[0].length)
    for (const s of seqs) {
        for (let i = 0; i < seqs[0].length; i++) {
            const idx = alphabet.indexOf(s[i])
            p[idx][i] = parseFloat((log(2, (f[idx][i] / 0.25))).toFixed(2).replace(/\.?0*$/, ''));
        }

    }

    return p
}
const getMotif = (mat, seqs, alphabet = "ACGT") => {
    let m = ""
    for (let j = 0; j < seqs[0].length; j++) {
        let max_c = mat[0][j]
        let max_ci = 0
        for (let i = 0; i < alphabet.length; i++) {
            if (mat[i][j] > max_c) {
                max_c = mat[i][j]
                max_ci = i
            }
        }
        m = m + alphabet[max_ci]

    }
    return m
}
const printMatInHtml = (mat) => {
    var result = "<table border='1' style='border-collapse:collapse'><tr>";
    result += "</tr>";
    var val;
    for (var i = 0; i < mat.length; i++) {
        result += "<tr>";
        for (var j = 0; j < mat[i].length; j++) {
            val = mat[i][j];
            result += "<td style='padding:1rem;'>" + val + "</td>";
        }
        result += "</tr>";
    }
    result += "</table>";

    return result;

}
//const seqs=[
//        "GCCGGAAGTG",
//        "ACCGGAAGCA",
//        "GCCGGATGTA",
//        "ACCGGAAGCT",
//        "ACCGGATATA",
//        "CCCGGAAGTG",
//        "ACAGGAAGTC",
//        "GCCGGATGCA",
//        "TCCGGAAGTA",
//        "ACAGGAAGCG",
//        "ACAGGATATG",
//        "TCCGGAAACC",
//        "ACAGGATATC",
//        "CAAGGACGAC",
//        "TCTGGACCCT",
//        ]


open.addEventListener('click', (e) => {
    e.preventDefault()
    dialog.setAttribute("open", "")
})

close.addEventListener('click', (e) => {
    e.preventDefault()
    dialog.removeAttribute("open")
})

const run = ()=>{
    let seqs = s.value.toUpperCase().trim("/s").split("\n")
    console.log(seqs)
    const c = matCount(seqs)
    const f = fqMat(seqs)
    const p = pwm(seqs)
    // const m = getMotif(seqs,p)


    console.log("matrice de contage")
    printMat(c)
    console.log("matrice de frequence")
    printMat(f)
    console.log("matrice de PWM")
    printMat(p)
    // console.log("Motif")
    // console.log(m)
    //

    const cmat = printMatInHtml(c)
    cm.innerHTML = cmat
    const fmat = printMatInHtml(f)
    fm.innerHTML = fmat
    const pmat = printMatInHtml(p)
    pm.innerHTML = pmat
}

btn.addEventListener('click', (e) => {
    e.preventDefault()
    run()

})

