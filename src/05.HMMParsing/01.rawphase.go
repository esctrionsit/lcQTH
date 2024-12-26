package main

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
)

func loadList(filepath string) (lines []string) {
	content, err := os.ReadFile(filepath)
	if err != nil {
		log.Fatal("[E] Error in loading file "+filepath+"\n", err)
	}
	for _, line := range strings.Split(string(content), "\n") {
		if len(line) > 0 {
			lines = append(lines, line)
		}
	}

	return
}

func loadCHRlen(filepath string) map[string]int {
	CHRlen := make(map[string]int)

	lines := loadList(filepath)
	for _, line := range lines {
		ele := strings.Split(line, "\t")
		lenInBp, _ := strconv.Atoi(ele[1])
		CHRlen[ele[0]] = lenInBp
	}

	return CHRlen
}

type jobPara struct {
	SM    string
	chrom string
}

type resPara struct {
	RawPhase      string
	ParentalRatio string
}

func jobPublisher(SMs []string, CHR string, jobs chan<- jobPara) {
	for _, SM := range SMs {
		newjob := jobPara{
			SM:    SM,
			chrom: CHR,
		}
		jobs <- newjob
	}
}

func phase(jobs <-chan jobPara, resultsChan chan<- resPara, idx int, CHRlen map[string]int, StepWidth int, WindowSize int) {
	for job := range jobs {
		// Load RR
		lines := loadList("../03.IdentifyRCR/ReadCovRegionList/" + job.SM + ".txt")
		rrstart := make([]int, 0, len(lines))
		for _, line := range lines {
			ele := strings.Split(line, "\t")
			if ele[0] == job.chrom {
				rrstartval, _ := strconv.Atoi(ele[1])
				rrstart = append(rrstart, rrstartval)
			}
		}

		// Load Rawphase
		lines = loadList("../04.HaplotypeRCR/Rawphase/" + job.SM + ".go.txt")
		var phaseseq string
		for _, line := range lines {
			ele := strings.Split(line, "\t")
			if ele[0] == job.chrom {
				phaseseq = ele[1]
			}
		}

		// clean
		lines = make([]string, 0, len(lines))

		// phase
		var PhasedSeq bytes.Buffer
		maxbins := int(math.Ceil(float64(CHRlen[job.chrom]) / (float64(StepWidth))))
		P1RatioLst := make([]string, 0, maxbins)
		P2RatioLst := make([]string, 0, maxbins)
		// IdenticalRatioLst := make([]string, 0, maxbins)
		UnknownRatioLst := make([]string, 0, maxbins)
		cache := []int{0, 0, 0, 0}
		cacheMap := map[string]int{"+": 0, ".": 1, "-": 2, "x": 3}
		RRcursor := 0

		for i := WindowSize; i <= CHRlen[job.chrom]; i += StepWidth {
			for RRcursor > 0 && rrstart[RRcursor-1] >= i-WindowSize {
				RRcursor--
			}
			for RRcursor < len(rrstart) && rrstart[RRcursor] < i-WindowSize {
				RRcursor++
			}

			for RRcursor < len(rrstart) && rrstart[RRcursor] < i && rrstart[RRcursor] >= i-WindowSize {
				cache[cacheMap[string(phaseseq[RRcursor])]] += 1
				RRcursor++
			}

			if cache[0]+cache[2] == 0 {
				PhasedSeq.WriteString(".")
				P1RatioLst = append(P1RatioLst, "0")
				P2RatioLst = append(P2RatioLst, "0")
				if cache[0]+cache[2]+cache[3] > 0 {
					UnknownRatioLst = append(UnknownRatioLst, fmt.Sprintf("%.2f", float64(cache[3])/float64(cache[0]+cache[2]+cache[3])*100))
				} else {
					UnknownRatioLst = append(UnknownRatioLst, "0")
				}
			} else {
				P1RatioLst = append(P1RatioLst, fmt.Sprintf("%.2f", float64(cache[0])/float64(cache[0]+cache[2]+cache[3])*100))
				P2RatioLst = append(P2RatioLst, fmt.Sprintf("%.2f", float64(cache[2])/float64(cache[0]+cache[2]+cache[3])*100))
				UnknownRatioLst = append(UnknownRatioLst, fmt.Sprintf("%.2f", float64(cache[3])/float64(cache[0]+cache[2]+cache[3])*100))

				PR := float64(cache[0]) / float64(cache[0]+cache[2]+cache[3]) * 100
				if PR <= 10 {
					PhasedSeq.WriteString("+")
				} else if PR > 10 && PR < 25 {
					PhasedSeq.WriteString("*")
				} else if PR >= 90 {
					PhasedSeq.WriteString("-")
				} else if PR > 75 && PR < 90 {
					PhasedSeq.WriteString("/")
				} else {
					PhasedSeq.WriteString("±")
				}
			}
			cache = []int{0, 0, 0, 0}
		}
		newresult := resPara{
			RawPhase:      job.SM + "\t" + PhasedSeq.String() + "\n",
			ParentalRatio: job.SM + "\t" + strings.Join(P1RatioLst, ",") + "\n" + job.SM + "\t" + strings.Join(P2RatioLst, ",") + "\n" + job.SM + "\t" + strings.Join(UnknownRatioLst, ",") + "\n",
		}
		resultsChan <- newresult
	}
}

func WriteFile(filePath string, content string) {
	// 路径，打开模式，文件权限
	file, err := os.OpenFile(filePath, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0666)
	if err != nil {
		fmt.Println("文件打开失败："+filePath, err)
	}
	defer file.Close()
	//写入文件时，使用带缓存的 *Writer
	write := bufio.NewWriter(file)
	write.WriteString(content)
	write.Flush()
}

func main() {
	args := os.Args
	if len(args) < 3 || args == nil {
		fmt.Println("Usage: " + args[0] + " <SMlst path> <CHRlen path> <MAXTHR>")
		os.Exit(1)
	}
	fmt.Println(args[1])
	Offsprings := loadList(args[1])
	CHRlen := loadCHRlen(args[2])
	// fmt.Println(Offsprings)
	MAXT, _ := strconv.Atoi(args[3])
	StepWidth, _ := strconv.Atoi(args[4])
	WindowSize, _ := strconv.Atoi(args[5])

	jobs := make(chan jobPara, 2*MAXT)
	resultsChan := make(chan resPara, 100)

	for w := 0; w < MAXT; w++ {
		// Run worker
		go phase(jobs, resultsChan, w, CHRlen, StepWidth, WindowSize)
	}

	for CHR, _ := range CHRlen {
		// fmt.Printf(CHR + "\n")

		go jobPublisher(Offsprings, CHR, jobs)

		var RawPhaseOut bytes.Buffer
		var ParentalRatioOut bytes.Buffer
		for r := 0; r < len(Offsprings); r++ {
			SMpattern := <-resultsChan
			RawPhaseOut.WriteString(SMpattern.RawPhase)
			ParentalRatioOut.WriteString(SMpattern.ParentalRatio)
		}

		WriteFile("data/"+CHR+".txt", RawPhaseOut.String())
		WriteFile("RawPhaseParentalRatio/"+CHR+".txt", ParentalRatioOut.String())
	}
	fmt.Println("\nFinished.")
}
