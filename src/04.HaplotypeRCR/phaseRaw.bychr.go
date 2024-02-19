package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io/ioutil"
	"log"
	"math/rand"
	"os"
	"os/exec"
	"strconv"
	"strings"
	"time"
)

var source = rand.NewSource(time.Now().UnixNano())

const charset = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

func RandString(length int) string {
	b := make([]byte, length)
	for i := range b {
		b[i] = charset[source.Int63()%int64(len(charset))]
	}
	return string(b)
}

func readReadRange(child string) (RR []string) {
	// open file
	f, err := os.Open("../03.IdentifyRCR/ReadCovRegionList/" + child + ".txt")
	if err != nil {
		log.Fatal(err)
	}
	// remember to close the file at the end of the program
	defer f.Close()

	// read the file line by line using scanner
	scanner := bufio.NewScanner(f)
	scanner.Scan()
	// i := 0
	for scanner.Scan() {
		// do something with a line
		// if i == 1000{
		// 	break
		// }
		RR = append(RR, strings.Replace(scanner.Text(), "\n", "", -1))
		// i++
	}
	return
}

type jobPara struct {
	SMlst      string
	chrom      string
	RRlst      []string
	querySites string
}

type resPara struct {
	chrom      string
	phaseTypes string
	// sideOut     string
}

func WriteQuerySites(filePath string, content string) {
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

func Pred(jobs <-chan jobPara, resultsChan chan<- resPara, idx int, vcfpath string, vcfSuffixes string, tmpPath string) {
	for job := range jobs {
		var Rtypes bytes.Buffer
		// var Sout bytes.Buffer

		randToken := RandString(5)
		WriteQuerySites(tmpPath + "/"+randToken, job.querySites)

		shell := "bcftools view " + vcfpath + "/" + job.chrom + vcfSuffixes + " -s " + job.SMlst + " -T " + tmpPath + "/" + randToken + " --threads 2 | sed 's/\\.\\/1/1\\/1/g'| sed 's/\\.\\/0/0\\/0/g' | bcftools query -f '%POS[\\t%GT:%GQ]\\n' -e 'N_MISSING>0 || MAF==0' "
		// fmt.Println(shell)
		cmd := exec.Command("/bin/bash", "-c", shell)
		stdout, err := cmd.StdoutPipe()
		if err != nil {
			log.Fatal(err)
		}
		if err := cmd.Start(); err != nil {
			log.Fatal(err)
		}
		data, err := ioutil.ReadAll(stdout)
		if err != nil {
			log.Fatal(err)
		}
		if err := cmd.Wait(); err != nil {
			log.Fatal(err)
		}
		stdout.Close()
		SNPs := make([]string, 0)
		if len(data) != 0 {
			SNPs = strings.Split(string(data), "\n")
		}
		data = nil
		cmd = exec.Command("rm", tmpPath + "/"+randToken)

		SNPidx := 0

		for RRidx := 0; RRidx < len(job.RRlst); RRidx++ {
			eleRR := strings.Split(string(job.RRlst[RRidx]), "\t")
			RRstart, _ := strconv.Atoi(eleRR[1]) // current read region start
			RRend, _ := strconv.Atoi(eleRR[2])

			CDC := []int{0, 0} // current diffcount P1, P2
			eleGT := make([]string, 0)
			SNPpos := -1

			for i := SNPidx; i >= 0; i-- {
				SNPidx = i
				eleGT = strings.Split(strings.Replace(SNPs[SNPidx], "\n", "", -1), "\t")
				if len(eleGT) < 4 {
					continue
				}
				SNPpos, _ = strconv.Atoi(eleGT[0])
				if SNPpos >= RRstart {
					continue
				} else {
					break
				}
			}

			for i := SNPidx; i < len(SNPs); i++ {
				SNPidx = i
				eleGT = strings.Split(strings.Replace(SNPs[SNPidx], "\n", "", -1), "\t")
				if len(eleGT) < 4 {
					continue
				}
				eleGQ := make([]float64, 4)
				for k := 1; k < 4; k++ {
					GQ := strings.Split(eleGT[k], ":")
					// if len(GQ) < 2 {
					// 	fmt.Println(GQ)
					// }
					eleGT[k] = GQ[0]
					GQval, err := strconv.ParseFloat(GQ[1], 64)
					if err != nil {
						eleGQ[k] = 0
					} else {
						eleGQ[k] = GQval
					}
				}
				SNPpos, _ = strconv.Atoi(eleGT[0])
				if SNPpos >= RRstart {
					if eleGT[1] == eleGT[2] || strings.Contains(eleGT[1], ".") || eleGT[1][0] != eleGT[1][2] || eleGQ[1] < 8 || strings.Contains(eleGT[2], ".") || eleGT[2][0] != eleGT[2][2] || eleGQ[2] < 8 || strings.Contains(eleGT[3], ".") || eleGT[3][0] != eleGT[3][2] || eleGQ[3] < 4 {
						// continue
					} else {
						if eleGT[3] != eleGT[1] {
							CDC[0]++
						}
						if eleGT[3] != eleGT[2] {
							CDC[1]++
						}
					}
					break
				}
			}

			if SNPidx < 0 || SNPpos < RRstart || SNPpos > RRend {
				Rtypes.WriteString(".")
				continue
			}

			for i := SNPidx; i < len(SNPs); i++ {
				SNPidx = i
				eleGT = strings.Split(strings.Replace(SNPs[SNPidx], "\n", "", -1), "\t")
				if len(eleGT) < 4 {
					continue
				}
				eleGQ := make([]float64, 4)
				for k := 1; k < 4; k++ {
					GQ := strings.Split(eleGT[k], ":")
					eleGT[k] = GQ[0]
					GQval, err := strconv.ParseFloat(GQ[1], 64)
					if err != nil {
						eleGQ[k] = 0
					} else {
						eleGQ[k] = GQval
					}
				}
				SNPpos, _ = strconv.Atoi(eleGT[0])
				if SNPpos <= RRend {
					if eleGT[1] == eleGT[2] || strings.Contains(eleGT[1], ".") || eleGT[1][0] != eleGT[1][2] || eleGQ[1] < 8 || strings.Contains(eleGT[2], ".") || eleGT[2][0] != eleGT[2][2] || eleGQ[2] < 8 || strings.Contains(eleGT[3], ".") || eleGT[3][0] != eleGT[3][2] || eleGQ[3] < 4 {
						// continue
					} else {
						if eleGT[3] != eleGT[1] {
							CDC[0]++
						}
						if eleGT[3] != eleGT[2] {
							CDC[1]++
						}
					}
				} else {
					break
				}
			}

			if CDC[0] > 2 && CDC[1] > 2 {
				Rtypes.WriteString("x")
			} else if CDC[0] == CDC[1] {
				Rtypes.WriteString(".")
			} else if CDC[0] <= 1 && CDC[1] <= 1 {
				Rtypes.WriteString(".")
			} else if CDC[0] < CDC[1] && CDC[0] <= 1 {
				Rtypes.WriteString("+")
			} else if CDC[0] > CDC[1] && CDC[1] <= 1 {
				Rtypes.WriteString("-")
			} else {
				Rtypes.WriteString(".")
			}
		}
		// fmt.Println(len(Rtypes.String()), len(job.RRlst))

		newjob := resPara{
			chrom:      job.chrom,
			phaseTypes: Rtypes.String(),
			// sideOut:  Sout.String(),
		}
		resultsChan <- newjob
	}
}

func jobPublisher(Parents []string, SM string, RRs []string, jobs chan<- jobPara, chrCount int) {
	Cchrom := ""
	Si := 0
	chromcount := 0
	var querySites bytes.Buffer
	for i, RR := range RRs {
		ele := strings.Split(string(RR), "\t")

		if Cchrom == "" {
			Cchrom = ele[0]
			chromcount += 1
		} else if Cchrom != ele[0] {
			chromcount += 1
			newjob := jobPara{
				SMlst:      Parents[0] + "," + Parents[1] + "," + SM,
				chrom:      Cchrom,
				RRlst:      RRs[Si:i],
				querySites: querySites.String(),
			}
			jobs <- newjob
			Cchrom = ele[0]
			Si = i
			querySites.Reset()
		}
		querySites.WriteString(ele[0] + "\t" + ele[1] + "\t" + ele[2] + "\n")
	}

	if chromcount != chrCount {
		log.Fatal("[E] " + SM + " does not have RR data for all " + strconv.Itoa(chrCount) + " chromosomes.")
	}

	newjob := jobPara{
		SMlst:      Parents[0] + "," + Parents[1] + "," + SM,
		chrom:      Cchrom,
		RRlst:      RRs[Si:],
		querySites: querySites.String(),
	}
	jobs <- newjob
}

func WriteFile(SM string, ResList []string, ChrList []string, folderPath string, backend string) {
	filePath := folderPath + "/" + SM + backend
	file, err := os.OpenFile(filePath, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0666)
	if err != nil {
		fmt.Println("文件打开失败", err)
	}
	defer file.Close()
	//写入文件时，使用带缓存的 *Writer
	write := bufio.NewWriter(file)
	for i, chrom := range ChrList {
		write.WriteString(chrom + "\t")
		write.WriteString(ResList[i] + "\n")
	}
	write.Flush()
}

func loadSMlst(filepath string) (SMlst []string) {
	// open file
	f, err := os.Open(filepath)
	if err != nil {
		log.Fatal(err)
	}
	// remember to close the file at the end of the program
	defer f.Close()

	// read the file line by line using scanner
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		// do something with a line
		if len(scanner.Text()) > 0 && scanner.Text() != "\n" {
			SMlst = append(SMlst, strings.Replace(scanner.Text(), "\n", "", -1))
		}
	}

	return
}

func main() {
	args := os.Args
	// if len(args) < 2 || args == nil {
	// 	fmt.Println("Usage: " + args[0] + " <SMlst path>")
	// 	os.Exit(1)
	// }
	Parents := []string{args[2], args[3]} // VCFid
	// fmt.Println(args[1])
	Offsprings := loadSMlst(args[1])
	// fmt.Println(Offsprings)
	MAXT, _ := strconv.Atoi(args[4])
	vcfPath := args[5]
	vcfSuffixes := args[6]
	tmpPath := args[7]
	CHRlen, _ := strconv.Atoi(args[8])

	jobs := make(chan jobPara, 2*MAXT)
	resultsChan := make(chan resPara, 100)

	for w := 0; w < MAXT; w++ {
		// Run worker
		go Pred(jobs, resultsChan, w, vcfPath, vcfSuffixes, tmpPath)
	}

	for _, SM := range Offsprings {
		// fmt.Println(SM)
		RRs := readReadRange(SM)
		go jobPublisher(Parents, SM, RRs, jobs, CHRlen)

		ResList := make([]string, MAXT)
		ChrList := make([]string, MAXT)
		// SOList := make([]string, 21)
		for r := 0; r < CHRlen; r++ {
			ResList[r] = "X"
			ChrList[r] = "X"
		}

		for r := 0; r < CHRlen; r++ {
			val := <-resultsChan
			ResList[r] = val.phaseTypes
			ChrList[r] = val.chrom
			// SOList[r] = val.sideOut
		}

		WriteFile(SM, ResList, ChrList, "Rawphase/", ".go.txt")
		// WriteFile(SM, SOList, ChrList, "MinCDC/", ".txt")
	}
	// fmt.Println("Finished.")
}
